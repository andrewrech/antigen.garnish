## ---- garnish_jaffa
#' Process gene fusions and return a data table for neoepitope prediction.
#'
#' Process [JAFFA](https://github.com/Oshlack/JAFFA) gene fusion `fasta` and `results.csv` output for neoepitope prediction using `garnish_affinity`.
#'
#' @param path Path to `jaffa_results.csv`.
#' @param db Character vector. One of `GRCm37`, `GRCm38`, `GRCh37`, or `GRCh38`.
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
#' @seealso \code{\link{garnish_affinity}}
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # load example jaffa output
#'    path <- "antigen.garnish_jaffa_results.csv" %T>%
#'      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.csv", .)
#'    fasta_path <- "antigen.garnish_jaffa_results.fasta" %T>%
#'      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)
#'
#'  # get predictions
#'    dt <- antigen.garnish::garnish_jaffa(path, db = "GRCm38", fasta_path) %>%
#'
#'  # add MHC info with antigen.garnish::list_mhc() compatible names
#'    .[, MHC := "H-2-Kb"] %>%
#'
#'  # get predictions
#'    antigen.garnish::garnish_affinity %>%
#'
#'  # summarize predictions
#'    antigen.garnish::garnish_summary %T>%
#'    print
#'}
#'
#' @export garnish_jaffa
#' @md

garnish_jaffa <- function(path, db, fasta_path){

  # check input
    if (missing(db) | !(db %chin% c("GRCm38",
                                    "GRCh38",
                                    "GRCm37",
                                    "GRCh37"))) stop("db must be one of GRCm37, GRCm38, GRCh37, or GRCh38)")
    if (!file.exists(fasta_path)) stop("fasta_path does not exist")
    if (!file.exists(path)) stop("path file does not exist")

    if (db == "GRCh38") host <- "aug2017.archive.ensembl.org"
    if (db == "GRCh37") host <- "grch37.ensembl.org"
    if (db == "GRCm38") host <- "feb2014.archive.ensembl.org"
    if (db == "GRCm37") host <- "may2012.archive.ensembl.org"

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


    dt[, fusion_uuid := parallel::mclapply(1:nrow(dt),
                                uuid::UUIDgenerate) %>% unlist]

  # split up gene fusion components

unfuse_genes <- function(col){

dtl <- parallel::mclapply(1:length(col), function(i){
               strsplit(col[i], split = ":", fixed = TRUE) %>% unlist
               })

gene_1 <- lapply(dtl, function(x){x[1]}) %>% unlist

gene_2 <- lapply(dtl, function(x){x[2]}) %>% unlist

      return(list(gene_1, gene_2))
     }

  dt[, c("gene_1", "gene_2") := unfuse_genes(`fusion genes`)]

  # incorporate fasta file input, split contigs at breakpoint
    fasta <- ShortRead::readFasta(fasta_path)

    seqs <- fasta@sread %>% data.table::as.data.table %>% .[, x %>% toupper] %>%
      data.table::as.data.table %>%
      data.table::setnames(".", "sread")

    contig <- fasta@id %>% data.table::as.data.table %>% .[, x] %>%
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


  # prep biomaRt
    if (grepl("GRCh", db)) bmds <- "hsapiens_gene_ensembl"
    if (grepl("GRCm", db)) bmds <- "mmusculus_gene_ensembl"

    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                             dataset = bmds,
                             host = host,
                            ensemblRedirect = FALSE)

    gn <- c(dt[, gene_1], dt[, gene_2]) %>% unique

  if (bmds == "hsapiens_gene_ensembl")
   query_n <- "hgnc_symbol"

  if (bmds == "mmusculus_gene_ensembl")
   query_n <- "mgi_symbol"

  # obtain transcript metadata
    var_dt <- biomaRt::getBM(
               attributes = c("ensembl_transcript_id",
                              query_n,
                              "ensembl_gene_id",
                              "description",
                              "chromosome_name",
                              "start_position",
                              "end_position",
                              "transcript_start",
                              "transcript_end",
                              "refseq_mrna"),
                             filters = query_n,
                             values = list(gn),
                             mart = mart) %>%
      data.table::as.data.table

 if (bmds == "mmusculus_gene_ensembl")
  var_dt %>% data.table::setnames("mgi_symbol", "external_gene_name")

  if (bmds == "hsapiens_gene_ensembl")
   var_dt %>% data.table::setnames("hgnc_symbol", "external_gene_name")

      trn <- var_dt[, ensembl_transcript_id %>% unique]

  # obtain coding and wt sequences

seqdtl <- lapply(c("coding", "peptide"), function(j){

    dt <- biomaRt::getSequence(type = "ensembl_transcript_id",
                               id = trn,
                               seqType = j,
                               mart = mart) %>%
                               data.table::as.data.table
                      return(dt)
                    })


    seqdt <- merge(seqdtl[[1]], seqdtl[[2]], by = "ensembl_transcript_id")
    var_dt <- merge(var_dt, seqdt, by = "ensembl_transcript_id")

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

tx_dt <- parallel::mclapply(tx_dt[, ensembl_gene_id %>% unique], function(i){
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

dt <- dt[, keep := parallel::mclapply(1:nrow(dt), function(i){
                  grepl(pattern = contig_1[i], x = coding_wt_1[i],
                        fixed = TRUE)
                  }) %>% unlist][keep == TRUE]

dt <- dt[, keep := parallel::mclapply(1:nrow(dt), function(i){
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
