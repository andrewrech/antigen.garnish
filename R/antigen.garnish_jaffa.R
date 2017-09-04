## ---- garnish_jaffa
#' Derive variants from JAFFA output. 
#'
#' Returns a data.table object of reformatted input from JAFFA results files for garnish_predictions.
#' Uses fasta and results.csv produced by JAFFA to generate all possible mutant fusion peptides predicted by
#' RNAseq reads for further downstream analysis by garnish_predictions.  
#'
#' @param path Full file path to jaffa_results.csv.
#' @param db Character vector. One of `GRCm37`, `GRCm38`, `GRCh37`, or `GRCh38`.
#' @param MHCdt Data.table object or full file path to a rio::import compatible file with format:
#' 
#'     Column name                 Example input
#'
#'     sample_id                   sample_1
#'     MHC                         HLA-A02:01 HLA-A03:01
#'                                 H-2-Kb H-2-Kb
#'                                 all
#'                                 HLA-DRB111:07 [second type]
#' @param fasta.file Full file path to jaffa_results.fasta.
#' @return A data table of mutant peptides including:
#' * **sample_id**: sample_id derived from jaffa_results.csv input.
#' * **pep_mut**: Full length mutant fusion peptides from all possible inframe exon-exon fusions predicted by JAFFA.
#' * **mutant_index**: Location of first mismatched AA between fusion peptide (pep_mut) and wt sequence of 5' gene of fusion.
#' * **fusion gene**: Colon separated gene names composing fusion derived from JAFFA output.
#' * **chrom1**: Chromosome of first fusion gene, format `chr[1-23,X,Y]`
#' * **base1**: From JAFFA output, base number on chromosome of fusion gene 1 breakpoint
#' * **chrom2**: Chromosome of second fusion gene, format `chr[1-23,X,Y]`
#' * **base2**: From JAFFA output, base number on chromosome of fusion gene 2 breakpoint
#' * **fusion_uuid**: Unique identifier for fusion event that will propogate through garnish_predictions.
#' * **pep_wt**: Full wt sequence of peptide from 5' fusion gene.
#' * **fus_tx**: Full sequence of combined coding transcripts of predicted fusion.
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'# run garnish_jaffa pipeline
#'
#'    # download example jaffa results files
#'    path <- "jaffa_results.csv"
#'    utils::download.file("http://get.rech.io/jaffa_results.csv", path)
#'    
#'    fasta.file <- "jaffa_results.fasta"
#'    utils::download.file("http://get.rech.io/jaffa_results.fasta", fasta.file)
#'    
#'    db <- "GRCm38"
#'    MHCdt <- data.table(sample_id = c("4662", "Abx7"),
#'                         MHC = "H-2-Kb H-2-Db H-2-IAb)
#'
#'    # extract variants
#'    dt <- antigen.garnish::garnish_jaffa(path, db, MHCdt, fasta.file) %>%
#'
#'    # predict neoepitopes
#'    antigen.garnish::garnish_predictions 
#'
#'    head(dt)
#'}
#'
#' @export garnish_jaffa
#' @md
  
garnish_jaffa <- function(path, db, MHCdt, fasta.file){
  if (missing(db) | !(db %chin% c("GRCm38",
"GRCh38", "GRCm37", "GRCh37"))) stop("Please provide argument db, character vector of length one
                        %chin% c(\"GRCh38\", \"GRCh37\", \"GRCm37\",\"GRCm38\")")
  
  if (missing(MHCdt)) stop("Please provide argument MHCdt, either a file or data.table object with 2 columns
                          1) sample_id, 2) MHC, a column of space separated MHC alleles for each unique sample,
                          to see supported alleles, list_mhc(), or use \"all\":
                          
                          sample_id                    MHC
                          A                            H-2-Kb H-2-IAb
                          B                            HLA-A*02:03
                          C                            all")
  if (missing(fasta.file) | !file.exists(fasta.file)) stop("Provide correct full path to jaffa_results.fasta")
  if (!((class(MHCdt) %chin% c("character", "data.table", "data.frame", "matrix")) %>% any)) stop("MHCdt must be a file path or data.table object.")
  
  if (class(MHCdt)[1] == "character"){
    if (!file.exits(MHCdt)) stop("MHCdt file not found.")
      
    MHCdt <- rio::import(MHCdt) %>% data.table::as.data.table
    
  }
  
  if(class(MHCdt)[1] == "data.frame"){
    
    MHCdt %<>% data.table::as.data.table
  }
  
  if (missing(path)) stop("Please provide a jaffa_results.csv")
  if (db == "GRCh38") host <- "ensembl.org"
  if (db == "GRCh37") host <- "grch37.ensembl.org"
  if (db == "GRCm38") host <- "ensembl.org"
  if (db == "GRCm37") host <- "may2012.archive.ensembl.org"
  if(!file.exists(path)) stop("Input must be full file path to a jaffa_results.csv")
  
  dt <- data.table::fread(path) %>% data.table::as.data.table
  
  if (!((c("sample", "fusion genes", "chrom1", "chrom2", "base1", "base2", "strand1", "strand2",
         "gap (kb)", "spanning pairs", "spanning reads", "inframe", "aligns", "rearrangement",
         "contig", "contig break", "classification", "known") %chin% names(dt)) %>% all)){
    stop("jaffa_results.csv may be malformed.  Check column names.")
  }

  dt <- dt %>% data.table::setnames("sample", "sample_id")

  dt <- dt[aligns == TRUE & rearrangement == TRUE & classification != "LowConfidence" & inframe == TRUE]
  
  dt[, fusion_uuid := parallel::mclapply(1:nrow(dt),
                                         uuid::UUIDgenerate) %>% unlist]
 
 ##split up gene fusion components 
  unfuse_genes <- function(col){
    
    dtl <- parallel::mclapply(1:length(col), function(i){
      strsplit(col[i], split = ":", fixed = TRUE) %>% unlist
    })
    
    gene_1 <- lapply(dtl, function(x){x[1]}) %>% unlist
    
    gene_2 <- lapply(dtl, function(x){x[2]}) %>% unlist
    
    return(list(gene_1, gene_2))
  }
  
  dt[, c("gene_1", "gene_2") := unfuse_genes(`fusion genes`)]
 
  ##incorporate fasta file input, split contigs at breakpoint 
  fasta <- ShortRead::readFasta(fasta.file)
  
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
  
  ##split up contigs into parts for each gene in fusion
  dt[, contig_1 := substr(x = sread, 1, `contig break`),
     by = 1:nrow(dt)] %>%
        .[, contig_2 := substr(x = sread, `contig break` + 1, nchar(sread)),
          by = 1:nrow(dt)]
  
  
  ##prep biomaRt 
  if (grepl("GRCh", db)) bmds <- "hsapiens_gene_ensembl"
  if (grepl("GRCm", db)) bmds <- "mmusculus_gene_ensembl"
  
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = bmds,
                           host = host)
  
  gn <- c(dt[, gene_1], dt[, gene_2]) %>% unique

    # obtain transcript metadata
    var_dt <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                            "external_gene_name", "ensembl_gene_id", "description", "chromosome_name",
                                            "start_position", "end_position", "transcript_start", "transcript_end",
                                            "transcript_length", "refseq_mrna"),
                             filters = "external_gene_name",
                             values = list(gn),
                             mart = mart) %>%
      data.table::as.data.table
  
  trn <- var_dt[, ensembl_transcript_id %>% unique]
    
  # obtain coding and wt sequences
  seqdtl <- lapply(c("coding", "peptide"), function(j){
    dt <- biomaRt::getSequence(type = "ensembl_transcript_id",
                         id = trn,
                         seqType = j,
                         mart = mart) %>% data.table::as.data.table
                    return(dt)
                  })
  
  
  seqdt <- merge(seqdtl[[1]], seqdtl[[2]], by = "ensembl_transcript_id")
  var_dt <- merge(var_dt, seqdt, by = "ensembl_transcript_id")
  
  ##add ensembl_gene_id to jaffa dt, then generate unique row for each transcript id
  ##for first gene of fusion and then for second
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
                data.table::setnames(c("external_gene_name", "ensembl_gene_id"),
                                     c("gene_2", "gene_id_2")),
              by = "gene_2")
  
  ##collapse various transcripts to then separate rows after
  tx_dt <- var_dt[, c("ensembl_gene_id", "ensembl_transcript_id")] %>% unique
  
  tx_dt <- parallel::mclapply(tx_dt[, ensembl_gene_id %>% unique], function(i){
    row <- tx_dt[ensembl_gene_id == i] %>%
          .[, txs := .[, ensembl_transcript_id] %>%
              paste0(collapse = ";")] %>%
                .[, c("ensembl_gene_id", "txs")]
    return(row)
        }) %>% data.table::rbindlist %>% unique
  
  ##merge dts with collapsed transcripts for genes 1 and 2
  
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
  
  ##separate/melt the dts for all combos of transcript id
  
  dt <- dt %>% tidyr::separate_rows("txs_1", sep = ";") %>%
            tidyr::separate_rows("txs_2", sep = ";")
  
  ##now merge with coding info from var_dt
  
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
  
  ##grep appropriate read/contig pattern and throw out transcript ids that could not have produced read
  
  dt <- dt[, keep := parallel::mclapply(1:nrow(dt), function(i){
                grepl(pattern = contig_1[i], x = coding_wt_1[i])
                }) %>% unlist][keep == TRUE] 
  dt <- dt[, keep := parallel::mclapply(1:nrow(dt), function(i){
                grepl(pattern = contig_2[i], x = coding_wt_2[i])
                  }) %>% unlist][keep == TRUE] 
  dt[, keep := NULL]
  
  ##make fusion transcripts
  dt[, fus1 := coding_wt_1 %>%
       stringr::str_extract(paste0("^.*", contig_1))] %>%
    .[, fus2 := coding_wt_2 %>%
        stringr::str_extract(paste0(contig_2, ".*$"))]
  dt[, fus_tx := paste0(fus1, fus2)]
  
  ##translate mutant peptides and get mutant index
  dt[, pep_fus := antigen.garnish::translate_cDNA(fus_tx) %>%
       stringr::str_replace(pattern = "\\*.*$", replacement = "")]
  ##fuzzy matches will return NA, must drop those rows to continue
  dt <- dt[!is.na(pep_fus)]
  
  ##get mutant indices
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
  
  ##clean up a little
  dt <- dt[, c("sample_id", "pep_fus", "mutant_index", "fusion genes", "chrom1", "base1",
               "chrom2", "base2", "fusion_uuid", "pep_wt_1", "fus_tx")] %>%
          data.table::setnames(c("pep_fus", "pep_wt_1"), c("pep_mut", "pep_wt"))
  
  dt <- merge(dt, MHCdt, by = "sample_id")
  
  return(dt)
}
