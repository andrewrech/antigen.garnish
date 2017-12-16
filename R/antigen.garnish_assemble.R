


## ---- list_mhc
#' List available MHC types using standard nomenclature
#'
#' @export list_mhc
#' @md

list_mhc <- function(){

system.file("extdata",
      "all_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      .$V1
  }



## ---- detect_mhc
#' Replace MHC with matching type in prediction tool format
#'
#' @param x Vector of HLA types named for program to convert to.
#' @param alleles Data table of 2 columns, 1. formatted allele names 2. prediction tool name (e.g. mhcflurry, mhcnuggets, netMHC).
#'
#' @export detect_mhc
#' @md

detect_mhc <- function(x, alleles){

  prog <- deparse(substitute(x))

    for (hla in (x %>% unique)){

      if (hla == "all") next

      # match hla allele alelle (end|(not longer allele))
      hla_re <- paste0(hla, "($|[^0-9])")

      allele <- alleles[type == prog, allele %>%
        .[stringi::stri_detect_regex(., hla_re) %>% which]] %>%
        # match to first in case of multiple matches
        # e.g. netMHCII
        .[1]

        if (allele %>% length == 0) allele <- NA

        x[x == hla] <- allele
        }

      return(x)
  }



## ---- get_metadata
#' Internal function to add metadata by `ensembl_transcript_id`
#'
#' @param dt Data table with `INFO` column.
#' @param humandb Character vector. One of `GRCh37` or `GRCh38`.
#' @param mousedb Character vector. One of `GRCm37` or `GRCm38`.
#'
#' @export get_metadata
#' @md

get_metadata <- function(dt,
                         humandb = "GRCh38",
                         mousedb = "GRCm38"){

  if (!"ensembl_transcript_id" %chin%
      (dt %>% names))
  stop("ensembl_transcript_id column missing")


  # set genome host
  if (!humandb %chin% c("GRCh37", "GRCh38")) stop("humandb set incorrectly")
  if (!mousedb %chin% c("GRCm37", "GRCm38")) stop("mousedb set incorrectly")
  if (humandb == "GRCh38") hhost <- "aug2017.archive.ensembl.org"
  if (humandb == "GRCh37") hhost <- "grch37.ensembl.org"
  if (mousedb == "GRCm38") mhost <- "ensembl.org"
  if (mousedb == "GRCm37") mhost <- "aug2017.archive.ensembl.org"

    # remove version suffix
    dt[, ensembl_transcript_id :=
      ensembl_transcript_id %>%
      stringr::str_replace("\\.[0-9]$", "")]

    bmds <- vector()

    if (dt[, ensembl_transcript_id %>%
        stringr::str_detect("ENSMUST")] %>%
          stats::na.omit %>%
          unique %>%
          any) bmds <- c(bmds, "mmusculus_gene_ensembl")

    if (dt[, ensembl_transcript_id %>%
        stringr::str_detect("ENST")] %>%
          stats::na.omit %>%
          unique %>%
          any) bmds <- c(bmds, "hsapiens_gene_ensembl")

    message("Obtaining cDNA and peptide sequences using biomaRt")

    var_dt <- lapply(bmds, function(i){

      if (i == "hsapiens_gene_ensembl") host <- hhost
      if (i == "mmusculus_gene_ensembl") host <- mhost

      mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                                 dataset = i,
                                                 host = host,
                                               ensemblRedirect = FALSE)

    if (i == "mmusculus_gene_ensembl"){
        trn <- dt[, ensembl_transcript_id %include% "ENSMUST" %>%
                  stats::na.omit %>%
                  unique]
              }

    if (i == "hsapiens_gene_ensembl"){
        trn <- dt[, ensembl_transcript_id %include% "ENST" %>%
                  stats::na.omit %>%
                  unique]
              }

   if (trn %>% length < 1) return(NULL)
   if (trn %>% length >= 1){

    # obtain transcript metadata

    # LPR - I have to take ref_seq mrna out because there are refseq mRNA that are deprecated with ensembl_transcript_ids
    # I have definitely had to remove it before.  Does it need to be in the package for anything specific?
    # Keeping it in causes the "allow cartesian" warning message to appear in the downstream merge with my 4662 VCFs.
      var_dt <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                    "external_gene_name", "ensembl_gene_id", "description", "chromosome_name",
                    "start_position", "end_position", "transcript_start", "transcript_end",
                    "transcript_length"),
                         filters = c("ensembl_transcript_id"),
                         values = list(trn),
                         mart = mart) %>%
            data.table::as.data.table

      # obtain transcript cDNA and peptide sequences

        seqdtl <- lapply(c("coding", "peptide"), function(j){
         biomaRt::getSequence(type = "ensembl_transcript_id",
                     id = trn,
                     seqType = j,
                     mart = mart) %>% data.table::as.data.table
      })

      seqdt <- merge(seqdtl[[1]], seqdtl[[2]], by = "ensembl_transcript_id")
      var_dt <- merge(var_dt, seqdt, by = "ensembl_transcript_id")
      }

    return(var_dt)

    }) %>% data.table::rbindlist
    dt <- merge(dt, var_dt, by = "ensembl_transcript_id")

    return(dt)
}



## ---- make_cDNA
#' Internal function to create cDNA from HGVS notation
#'
#' @param dt Data table with INFO column.
#'
#' @export make_cDNA
#' @md

make_cDNA <- function(dt){

  if (!c("cDNA_type",
         "coding",
         "cDNA_locs",
         "cDNA_locl",
         "cDNA_seq"
         ) %chin%
      (dt %>% names) %>% all)
  stop("dt is missing columns")

  # initialize and set types

    dt[, coding_mut := coding]

    for (i in c("cDNA_locs", "cDNA_locl")){
      set(dt, j = i, value = dt[, get(i) %>%
          as.integer])
    }

  # remove if nonsensical loc
    dt %<>% .[, coding_nchar := coding %>% nchar]

    dt %<>%
     .[!cDNA_locs > coding_nchar &
       !cDNA_locl > coding_nchar]

    dt[, coding_nchar := NULL]



  # paste substr around changed bases
  # substr is vectorized C, fast
  # mut sub nmers that are really wt
    # will be filtered after subpep generation


    # handle single base changes

      dt[cDNA_type == ">",
      coding_mut := coding_mut %>%
        { paste0(
           substr(., 1, cDNA_locs-1),
           cDNA_seq,
           substr(.,
                  cDNA_locs + 1,
                  nchar(.))) } ]

    # handle deletions
      # indices are inclusive
      # regexpr match for del

      dt[cDNA_type %like% "del",
      coding_mut := coding_mut %>%
        { paste0(
           substr(., 1, cDNA_locs-1),
           substr(.,
                  cDNA_locl + 1,
                  nchar(.))) } ]

    # handle insertions
      # regexpr match for del

      dt[cDNA_type %like% "ins",
      coding_mut := coding_mut %>%
        { paste0(
           substr(., 1, cDNA_locs),
           cDNA_seq,
           substr(.,
                  cDNA_locs + 1,
                  nchar(.))) } ]

      }



## ---- get_snpeff
#' Internal function to extract SnpEff annotation information to a data table.
#'
#' @param dt Data table with INFO column from a SnpEff-annotated VCF file.
#' @export get_snpeff
#' @md

get_snpeff <- function(dt){

    if (!"INFO" %chin% (dt %>% names)) stop("dt must contain INFO column")

    dt[, se := INFO %>%
      stringr::str_extract("ANN.*") %>%
      stringr::str_replace("ANN=[^\\|]+\\|", "")]

    # add a variant identifier
    suppressWarnings(dt[, snpeff_uuid :=
                  lapply(1:nrow(dt),
                  uuid::UUIDgenerate) %>% unlist])

    # abort if no variants passed filtering
    if (dt %>% nrow < 1) return(NULL)

    # spread SnpEff annotation over rows
    dt %>% tidyr::separate_rows("se", sep = ",")

    # extract info from snpeff annotation
      dt[, effect_type := se %>%
          stringr::str_extract("^[a-z0-9][^\\|]+")]
      dt[, ensembl_transcript_id := se %>%
          stringr::str_extract("(?<=\\|)(ENSMUST|ENST)[0-9]+")]
      dt[, ensembl_gene_id := se %>%
          stringr::str_extract("(?<=\\|)(ENSMUSG|ENSG)[0-9.]+(?=\\|)")]
      dt[, protein_change := se %>%
          stringr::str_extract("p\\.[^\\|]+")]
      dt[, cDNA_change := se %>%
          stringr::str_extract("c\\.[^\\|]+")]
      dt[, protein_coding := se %>%
          stringr::str_detect("protein_coding")]

      return(dt)

        }



## ---- extract_cDNA
#' Internal function to extract cDNA changes from HGVS notation
#'
#' @param dt Data table with INFO column.
#'
#' @export extract_cDNA
#' @md

extract_cDNA <- function(dt){

  # check required cols

  if (!c("cDNA_change") %chin%
        (dt %>% names) %>% all)
    stop("cDNA_change missing")


  # extract HGVS DNA nomenclature
      # dups are just ins
      dt[, cDNA_change := cDNA_change %>%
            stringr::str_replace_all("dup", "ins")]

      dt[, cDNA_locs := cDNA_change %>%
            stringr::str_extract("[0-9]+") %>%
            as.integer]
      dt[, cDNA_locl := cDNA_change %>%
            stringr::str_extract("(?<=_)[0-9]+") %>%
            as.integer]
        # make cDNA_locl cDNA_locs for single bases
             dt[cDNA_locl %>% is.na, cDNA_locl := cDNA_locs]

      dt[, cDNA_type := cDNA_change %>%
            stringr::str_extract_all("[a-z]{3}|>") %>%
            unlist %>%
            paste(collapse = ""),
            by = 1:nrow(dt)]
      dt[, cDNA_seq := cDNA_change %>%
            stringr::str_extract("[A-Z]+$")]

  # filter out extraction errors
  # https://github.com/andrewrech/antigen.garnish/issues/3
    dt %<>% .[!cDNA_locs %>% is.na &
       !cDNA_locl %>% is.na &
       !cDNA_type %>% is.na]

  }



## ---- translate_cDNA
#' Internal function to translate cDNA to peptides
#'
#' @param v cDNA character vector without ambiguous bases.
#'
#' @export translate_cDNA
#' @md

translate_cDNA <- function(v){

parallel::mclapply(v, function(p){

     # protect vector length
       tryCatch({

          p %>%
            Biostrings::DNAString() %>%
            Biostrings::translate() %>%
                as.vector %>%
                as.character %>%
                paste(collapse = "")

}, error = function(e){
              return(NA)
              })

        }) %>% unlist
  }
