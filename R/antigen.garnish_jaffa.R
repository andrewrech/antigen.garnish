
#' Process gene fusions and return a data table for neoantigen prediction.
#'
#' Process [JAFFA](https://github.com/Oshlack/JAFFA) gene fusion `fasta` and `results.csv` output for neoantigen prediction using `garnish_affinity`.
#'
#' @param path Path to `jaffa_results.csv`.
#' @param db One of "GRCm38" or "GRCh38", murine and human reference genomes respectively.
#' @param fasta_path Path to `jaffa_results.fasta`.
#' @return A data table of mutant peptides, including:
#' * **sample_id**: sample id
#' * **pep_mut**: mutant fusion peptide
#' * **mutant_index**: index of mutant peptide
#' * **fusion gene**: colon-separated gene names for the fusion product
#' * **chrom1**: chromosome of first fusion gene
#' * **base1**: base number of first fusion gene breakpoint
#' * **chrom2**: chromosome of second fusion gene
#' * **base2**: base number of second fusion gene breakpoint
#' * **fusion_uuid**: unique identifier
#' * **pep_wt**: wt cDNA sequence of peptide from 5' fusion gene.
#' * **fus_tx**: cDNA sequence of predicted fusion product
#'
#' @seealso \code{\link{garnish_variants}}
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_summary}}
#'
#' @examples
#'\dontrun{
#'
#'  # load example jaffa output
#'  dir <- system.file(package = "antigen.garnish") %>%
#'    file.path(., "extdata/testdata")
#'
#'    path <- "antigen.garnish_jaffa_results.csv" %>%
#'      file.path(dir, .)
#'    fasta_path <- "antigen.garnish_jaffa_results.fasta" %>%
#'      file.path(dir, .)
#'
#'  # get predictions
#'    dt <- garnish_jaffa(path, db = "GRCm38", fasta_path) %>%
#'
#'  # add MHC info with list_mhc() compatible names
#'    .[, MHC := "H-2-Kb"] %>%
#'
#'  # get predictions
#'    garnish_affinity %>%
#'
#'  # summarize predictions
#'    garnish_summary %T>%
#'    print
#'}
#'
#' @export garnish_jaffa
#' @md

garnish_jaffa <- function(path, db, fasta_path){

  # check input
    if (!file.exists(fasta_path)) stop("fasta_path does not exist")
    if (!file.exists(path)) stop("path file does not exist")
    if (!db %chin% c("GRCm38", "GRCh38")) stop("db must be either GRCm38 (murine) or GRCh38 (human).")

    # magrittr version check, this will not hide the error, only the NULL return on successful exit
    invisible(check_dep_versions())

    dt <- data.table::fread(path) %>%
      data.table::as.data.table

    if (!((c("sample",
             "fusion genes",
             "chrom1",
             "chrom2",
             "base1",
             "base2",
             "strand1",
             "strand2",
             "gap (kb)",
             "spanning pairs",
             "spanning reads",
             "inframe",
             "aligns",
             "rearrangement",
             "contig",
             "contig break",
             "classification",
             "known") %chin% names(dt)) %>% all)){
      stop("jaffa_results.csv has incorrect column names.")
    }

    dt <- dt %>% data.table::setnames("sample", "sample_id")

    dt <- dt[aligns == TRUE &
             rearrangement == TRUE &
             classification != "LowConfidence" &
             inframe == TRUE]


    dt[, fusion_uuid := uuid::UUIDgenerate(), by = 1:nrow(dt)]

  # split up gene fusion components

unfuse_genes <- function(col){

dtl <- lapply(1:length(col), function(i){
               strsplit(col[i], split = ":", fixed = TRUE) %>% unlist
               })

gene_1 <- lapply(dtl, function(x){x[1]}) %>% unlist

gene_2 <- lapply(dtl, function(x){x[2]}) %>% unlist

      return(list(gene_1, gene_2))
     }

  dt[, c("gene_1", "gene_2") := unfuse_genes(`fusion genes`)]

  # incorporate fasta file input, split contigs at breakpoint
    fasta <- Biostrings::readDNAStringSet(fasta_path)

    seqs <- fasta %>% as.character %>%
      data.table::as.data.table %>%
      data.table::setnames(".", "sread")

    contig <- fasta %>% as.character %>% names %>%
              data.table::as.data.table %>%
      data.table::setnames(".", "contig")

    contig_dt <- data.table::data.table(seqs, contig) %>%
        .[, contig := contig %>% stringr::str_extract("(?<=(\\-{3})).*$") %>%
            stringr::str_extract("(?<=(\\-{3})).*$")]


    dt <- merge(dt,
          contig_dt,
          by = "contig")

  # split up contigs into parts for each gene in fusion
    dt[, contig_1 := substr(x = sread, 1, `contig break`),
       by = 1:nrow(dt)] %>%
          .[, contig_2 := substr(x = sread, `contig break` + 1, nchar(sread)),
            by = 1:nrow(dt)]


  # detect/set AG_DATA_DIR environmental variable
  check_pred_tools()

  # load metadata
   metafile <- file.path(Sys.getenv("AG_DATA_DIR"), "/", "GRChm38_meta.RDS")

  if (!file.exists(metafile)){

    ag_data_err()

  }

    var_dt <- readRDS(metafile)

    # get species specific table here.
    if (db == "GRCh38") var_dt <- var_dt[ensembl_gene_id %likef% "ENSG"]
    if (db == "GRCm38") var_dt <- var_dt[ensembl_gene_id %likef% "ENSMUSG"]

    # toss what we don't need here
    var_dt <- var_dt[external_gene_name %chin%
    c(dt[, gene_1 %>% unique], dt[, gene_2 %>% unique])]
  # add ensembl_gene_id to jaffa dt
  # generate unique row for each transcript id
  # for first gene of fusion and then for second
    dt <- merge(dt,
              var_dt[external_gene_name %chin% dt[, gene_1 %>% unique],
                     c("external_gene_name", "ensembl_gene_id")] %>%
                       unique %>%
              data.table::setnames(c("external_gene_name", "ensembl_gene_id"),
                                   c("gene_1", "gene_id_1")),
                by = "gene_1")

    dt <- merge(dt,
              var_dt[external_gene_name %chin% dt[, gene_2 %>% unique],
                     c("external_gene_name", "ensembl_gene_id")] %>%
                       unique %>%
             data.table::setnames(c("external_gene_name",
                                    "ensembl_gene_id"),
                                  c("gene_2",
                                    "gene_id_2")),
                 by = "gene_2")

  # collapse transcripts to then separate rows after
    tx_dt <- var_dt[, c("ensembl_gene_id", "ensembl_transcript_id")] %>% unique

tx_dt <- lapply(tx_dt[, ensembl_gene_id %>% unique], function(i){
      row <- tx_dt[ensembl_gene_id == i] %>%
            .[, txs := .[, ensembl_transcript_id] %>%
                paste0(collapse = ";")] %>%
                  .[, c("ensembl_gene_id", "txs")]
      return(row)
          }) %>% data.table::rbindlist %>% unique

  # merge dts with collapsed transcripts for genes 1 and 2
    dt <- merge(dt,
                tx_dt[ensembl_gene_id %chin% dt[, gene_id_1 %>% unique]] %>%
                  data.table::setnames(c("ensembl_gene_id", "txs"),
                                       c("gene_id_1", "txs_1")),
                  by = "gene_id_1")

    dt <- merge(dt,
                tx_dt[ensembl_gene_id %chin% dt[, gene_id_2 %>% unique]] %>%
                  data.table::setnames(c("ensembl_gene_id", "txs"),
                                       c("gene_id_2", "txs_2")),
                by = "gene_id_2")

  # separate/melt the dts for all combinations of transcript id
    dt <- dt %>% tidyr::separate_rows("txs_1", sep = ";") %>%
                 tidyr::separate_rows("txs_2", sep = ";")

  # now merge with coding info from var_dt
    dt <- merge(dt,
                var_dt[ensembl_gene_id %chin% dt[, gene_id_1 %>% unique],
                       c("ensembl_transcript_id", "coding", "peptide")] %>%
                  unique %>%
                  data.table::setnames(c("ensembl_transcript_id", "coding", "peptide"),
                                       c("txs_1", "coding_wt_1", "pep_wt_1")),
                 by = "txs_1")

     dt <- merge(dt,
              var_dt[ensembl_gene_id %chin% dt[, gene_id_2 %>% unique],
                     c("ensembl_transcript_id", "coding", "peptide")] %>%
                unique %>%
                data.table::setnames(c("ensembl_transcript_id", "coding", "peptide"),
                                     c("txs_2", "coding_wt_2", "pep_wt_2")),
              by = "txs_2")

  # grep appropriate read/contig pattern
  # remove transcript ids that could not have produced read

dt <- dt[, keep := lapply(1:nrow(dt), function(i){
                  grepl(pattern = contig_1[i], x = coding_wt_1[i],
                        fixed = TRUE)
                  }) %>% unlist][keep == TRUE]

dt <- dt[, keep := lapply(1:nrow(dt), function(i){
                  grepl(pattern = contig_2[i], x = coding_wt_2[i],
                        fixed = TRUE)
                    }) %>% unlist][keep == TRUE]
    dt[, keep := NULL]

  # make fusion transcripts
    dt[, fus1 := coding_wt_1 %>%
         stringr::str_extract(paste0("^.*", contig_1))] %>%
      .[, fus2 := coding_wt_2 %>%
          stringr::str_extract(paste0(contig_2, ".*$"))]
    dt[, fus_tx := paste0(fus1, fus2)]

  # translate mutant peptides
  # get mutant index
    dt[, pep_fus := translate_cDNA(fus_tx) %>%
         stringr::str_replace(pattern = "\\*.*$", replacement = "")]

  # fuzzy matches will return NA, must drop these
    dt <- dt[!is.na(pep_fus)]

  # get mutant indices
    suppressWarnings({
    dt[, mutant_index :=
    {
      (pep_wt_1 %>%
         strsplit(split = "", fixed = TRUE) %>%
         unlist) !=
        (pep_fus %>%
           strsplit(split = "", fixed = TRUE) %>%
           unlist)
    } %>%
      which %>% .[1], by = 1:nrow(dt)]
      })

  # clean up columns
    dt <- dt[, c("sample_id",
                 "pep_fus",
                 "mutant_index",
                 "fusion genes",
                 "fusion_uuid",
                 "chrom1",
                 "base1",
                 "chrom2",
                 "base2",
                 "fus_tx",
                 "pep_wt_1")] %>%
            data.table::setnames(c("pep_fus", "pep_wt_1"),
                                 c("pep_mut", "pep_gene_1"))

  return(dt)
}
