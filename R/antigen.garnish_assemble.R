#' Return a data table of available MHC types for all prediction tools.
#'
#' @export list_mhc
#' @md

list_mhc <- function() {
  dt <- system.file("extdata",
    "all_alleles.txt",
    package = "antigen.garnish"
  ) %>%
    data.table::fread(header = FALSE, sep = "\t") %>%
    data.table::setnames("V1", "MHC") %>%
    .[, species := "human"] %>%
    .[MHC %like% "H-2", species := "mouse"] %>%
    .[, class := "II"] %>%
    .[MHC %like% "(H-2-[A-Z][a-z])|(HLA-[ABCE]\\*)", class := "I"]

  return(dt)
}


#' Internal function to replace MHC with matching prediction tool MHC syntax
#'
#' @param x Vector of HLA types named for program to convert to.
#' @param alleles Data table of 2 columns, 1. formatted allele names 2. prediction tool name (e.g. mhcflurry, netMHC).
#'
#' @noRd

detect_mhc <- function(x, alleles) {
  prog <- deparse(substitute(x))

  for (hla in (x %>% unique())) {
    if (hla %like% "all") next

    # match hla allele alelle (end|(not longer allele))
    hla_re <- paste0(hla, "($|[^0-9])")

    allele <- alleles[type == prog, allele %>%
      .[stringr::str_detect(., hla_re) %>% which()]] %>%
      # match to first in case of multiple matches
      # e.g. netMHCII
      .[1]

    if (allele %>% length() == 0) allele <- NA

    x[x == hla] <- allele
  }

  return(x)
}


#' Internal function to add metadata by `transcript_id`
#'
#' @param dt Data table with `INFO` column.
#'
#' @noRd

get_metadata <- function(dt) {
  if (!"transcript_id" %chin%
    (dt %>% names())) {
    stop("transcript_id column missing")
  }

  message("Reading local transcript metadata.")

  # detect/set AG_DATA_DIR environmental variable
  check_pred_tools()

  # load metadata
  metafile <- file.path(Sys.getenv("AG_DATA_DIR"), "/GRChm38_meta.RDS")

  if (!file.exists(metafile)) {
    .ag_data_err()
  }

  var_dt <- readRDS(metafile)

  dt <- merge(dt, var_dt, by = "transcript_id")

  rm(var_dt)

  return(dt)
}


#' Internal function to create cDNA from HGVS notation
#'
#' @param dt Data table with INFO column.
#'
#' @noRd

make_cDNA <- function(dt) {
  if (!c(
    "cDNA_type",
    "coding",
    "cDNA_locs",
    "cDNA_locl",
    "cDNA_seq"
  ) %chin%
    (dt %>% names()) %>% all()) {
    stop("dt is missing columns")
  }

  # initialize and set types

  dt[, coding_mut := coding]

  for (i in c("cDNA_locs", "cDNA_locl")) {
    set(dt, j = i, value = dt[, get(i) %>%
      as.integer()])
  }

  # remove if nonsensical loc
  dt %<>% .[, coding_nchar := coding %>% nchar()]

  dt %<>%
    .[!cDNA_locs > coding_nchar &
      !cDNA_locl > coding_nchar]

  dt[, coding_nchar := NULL]

  # paste substr around changed bases
  # substr is vectorized C, fast
  # mut sub nmers that are really wt
  # will be filtered after subpep generation


  # handle single base changes

  dt[
    cDNA_type == ">",
    coding_mut := coding_mut %>% {
      paste0(
        substr(., 1, cDNA_locs - 1),
        cDNA_seq,
        substr(
          .,
          cDNA_locs + 1,
          nchar(.)
        )
      )
    }
  ]

  # handle deletions
  # indices are inclusive
  # regexpr match for del

  dt[
    cDNA_type %like% "del",
    coding_mut := coding_mut %>% {
      paste0(
        substr(., 1, cDNA_locs - 1),
        substr(
          .,
          cDNA_locl + 1,
          nchar(.)
        )
      )
    }
  ]

  # handle insertions
  # regexpr match for del
  # this is off register by one for delins
  dt[
    cDNA_type == "ins",
    coding_mut := coding_mut %>% {
      paste0(
        substr(., 1, cDNA_locs),
        cDNA_seq,
        substr(
          .,
          cDNA_locs + 1,
          nchar(.)
        )
      )
    }
  ]

  # for delins, the cDNA_locs is a deleted base, so true s is now 1 off
  dt[
    cDNA_type == "delins",
    coding_mut := coding_mut %>% {
      paste0(
        substr(., 1, cDNA_locs - 1),
        cDNA_seq,
        substr(
          .,
          cDNA_locs,
          nchar(.)
        )
      )
    }
  ]
}


#' Internal function to extract a data table of variants with `INFO` fields in columns.
#'
#' @param vcf vcfR object to extract data from.
#'
#' @return Data table of variants with `INFO` fields in columns.
#'
#' @noRd

get_vcf_info_dt <- function(vcf) {
  if (vcf %>% class() %>% .[1] != "vcfR") {
    stop("vcf input is not a vcfR object.")
  }

  dt <- vcf@fix %>% data.table::as.data.table()

  if (!"INFO" %chin% (dt %>% names())) {
    stop("Error parsing input file INFO field.")
  }

  # loop over INFO field
  # tolerant of variable length and content
  v <- dt[, INFO %>% paste(collapse = ";@@@=@@@;")]
  vd <- v %>%
    stringr::str_replace_all("(?<=;)[^=;]+", "") %>%
    stringr::str_replace_all(stringr::fixed(";=@"), ";@") %>%
    stringr::str_replace_all(stringr::fixed("@;="), "@;")

  vn <- v %>%
    stringr::str_replace_all("(?<==)[^;]+", "") %>%
    stringr::str_replace_all(stringr::fixed("="), "")

  vd %<>% strsplit("@@@")
  vn %<>% strsplit(("@@@"))

  idt <- lapply(1:length(vn[[1]]), function(i) {
    v <- vd[[1]][i] %>%
      strsplit(";") %>%
      .[[1]]
    names(v) <- vn[[1]][i] %>%
      strsplit(";") %>%
      .[[1]]

    return(v %>% as.list())
  }) %>%
    rbindlist(fill = TRUE, use.names = TRUE) %>%
    .[, V1 := NULL] # remove extra final column created during parsing

  # remove extra colum

  if (
    (dt %>% nrow()) !=
      (idt %>% nrow())
  ) {
    stop("Error parsing input file INFO field.")
  }

  dt <- cbind(dt, idt)

  return(dt)
}


#' Internal function to extract a data table from vcfR `vcf` object fields.
#'
#' @param vcf vcfR object to extract data from.
#'
#' @return Data table of variants with sample level fields in columns.
#'
#' @noRd

get_vcf_sample_dt <- function(vcf) {
  if (vcf %>% class() %>% .[1] != "vcfR") {
    stop("vcf input is not a vcfR object.")
  }

  dt <- vcf@gt %>% data.table::as.data.table()

  if (!"FORMAT" %chin% (dt %>% names())) {
    stop("Error parsing input file sample level info.")
  }

  names <- vcf@gt %>%
    attributes() %>%
    .$dimnames %>%
    unlist() %exclude%
    "FORMAT"

  # loop over sample level data
  # tolerant of variable length and content

  # This used to  be doubly parallelized over samples in the vcf
  # would get weird recycling consistently with test data so I pulled this off.
  # this is will not substantially affect speed, only 2 samples to fork anyways.
  idt <- lapply(names, function(n) {
    ld <- dt[, get(n)]
    ln <- dt[, FORMAT]

    idt <- lapply(1:length(ln), function(i) {

      #  account for format columns longer than provided values, as in NORMAL vs. TUMOR with mutant
      # only one value is provided in tumor column for some values, like SA_POST_PROB  etc.
      nl <- ln[i] %>%
        strsplit(":") %>%
        .[[1]]
      l <- ld[i] %>%
        strsplit(":") %>%
        .[[1]]

      names(l) <- nl[1:length(l)]

      return(l %>% as.list())
    }) %>% rbindlist(fill = TRUE, use.names = TRUE)

    idt %>% data.table::setnames(
      idt %>% names(),
      idt %>% names() %>% paste0(n, "_", .)
    )
    return(idt)
  }) %>% do.call(cbind, .)

  # assign ref and alt to individual columns

  for (i in (idt %>% names() %include% "_AD$")) {
    colNum <- idt[, get(i) %>%
      data.table::tstrsplit(",") %>% length()]

    if (colNum == 2) {
      idt[, paste0(i, c("_ref", "_alt")) := get(i) %>%
        data.table::tstrsplit(",")]
    }
    if (colNum == 3) {
      idt[, paste0(i, c("_ref", "_alt", "_alt2")) := get(i) %>%
        data.table::tstrsplit(",")]
    }
    if (colNum == 4) {
      idt[, paste0(i, c("_ref", "_alt", "_alt2", "_alt3")) := get(i) %>%
        data.table::tstrsplit(",")]
    }
    if (colNum == 5) {
      idt[, paste0(i, c("_ref", "_alt", "_alt2", "_alt3", "_alt4")) := get(i) %>%
        data.table::tstrsplit(",")]
    }
    if (colNum > 5) {
      stop("Parsing VCFs with greater than 4 alternative alleles is not supported.")
    }

    set(idt, j = i, value = NULL)
  }

  if (
    (dt %>% nrow()) !=
      (idt %>% nrow())
  ) {
    stop("Error parsing input file INFO field.")
  }

  dt <- cbind(dt, idt)

  return(dt)
}


#' Internal function to extract SnpEff annotation information to a data table.
#'
#' @param dt Data table with character vector `ANN` column, from `vcf` file.
#'
#' @return Data table with the `ANN` column parsed into additional rows.
#'
#' @noRd

get_vcf_snpeff_dt <- function(dt) {
  if (!"ANN" %chin% (dt %>% names())) {
    stop("Error parsing input file ANN field from SnpEff.")
  }

  # add a variant identifier

  len <- nrow(dt)
  suppressWarnings(dt[, snpeff_uuid := uuid::UUIDgenerate(use.time = FALSE, n = len)])

  # spread SnpEff annotation over rows
  # transform to data frame intermediate to avoid
  # data.table invalid .internal.selfref
  dt %<>%
    # first, remove ALT allele annotation to standardize form across SnpEff annotations in each row
    .[, ANN := ANN %>% stringr::str_replace("^(ANN)?(=)?[ACTNG]*\\|", "")] %>%
    # second, separate SnpEff annotations into one per row
    tidyr::separate_rows("ANN", sep = ",([ACTG]+)?\\|") %>%
    data.table::as.data.table(.)

  dt[, transcript_id := ANN %>%
    stringr::str_extract("(?<=\\|)(ENSMUST|ENST)[0-9]+(\\.[0-9]+)?")]

  # fall back to RefSeq
  dt[transcript_id %>% is.na(), transcript_id := ANN %>%
    stringr::str_extract("(?<=\\|)NM_[0-9]+\\.[0-9]+")]

  # the nucleotide header must be removed from the ANN field to parse correctly
  # because some ANN fields do not have a header, resulting in the wrong
  # number of columns

  annSpecNames <- c(
    "effect_type",
    "putative_impact",
    "gene",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_bioptype",
    "exon_intron_rank",
    "cDNA_change",
    "protein_change",
    "cDNA_position_cDNA_len",
    "CDS_position_CDS_len",
    "Protein_position_Protein_len",
    "Distance_to_feature"
  )

  snpEff <- dt[, ANN %>% data.table::tstrsplit("\\|")]

  # also parse error column if it exists
  if (snpEff %>% length() == 15) {
    annSpecNames %<>% c("ERRORS_WARNINGS_INFO")
  }

  snpEff %>% data.table::setnames(annSpecNames)

  dt <- try(cbind(dt, snpEff))

  if (any(class(dt) %like% "error")) {
    print(str(dt))
    stop("error parsing VCF SnpEff annotations")
  }
  dt[, protein_coding := FALSE]
  dt[protein_change != "", protein_coding := TRUE]

  return(dt)
}


#' Internal function to extract cDNA changes from HGVS notation
#'
#' @param dt Data table with INFO column.
#'
#' @noRd

extract_cDNA <- function(dt) {

  # check required cols

  if (!c("cDNA_change") %chin%
    (dt %>% names()) %>% all()) {
    stop("cDNA_change missing")
  }


  # extract HGVS DNA nomenclature
  # dups are just ins
  dt[, cDNA_change := cDNA_change %>%
    stringr::str_replace_all("dup", "ins")]

  dt[, cDNA_locs := cDNA_change %>%
    stringr::str_extract("[0-9]+") %>%
    as.integer()]
  dt[, cDNA_locl := cDNA_change %>%
    stringr::str_extract("(?<=_)[0-9]+") %>%
    as.integer()]
  # make cDNA_locl cDNA_locs for single bases
  dt[cDNA_locl %>% is.na(), cDNA_locl := cDNA_locs]

  dt[, cDNA_type := cDNA_change %>%
    # make cDNA deletes then inserts, if we have indel variants it will work
    # the pattern pulled will be delins
    stringr::str_extract_all("[a-z]{3}|>") %>%
    unlist() %>%
    paste(collapse = ""),
  by = 1:nrow(dt)
  ]
  dt[, cDNA_seq := cDNA_change %>%
    stringr::str_extract("[A-Z]+$")]

  # filter out extraction errors
  # https://github.com/andrewrech/antigen.garnish/issues/3
  dt %<>% .[!cDNA_locs %>% is.na() &
    !cDNA_locl %>% is.na() &
    !cDNA_type %>% is.na()]
}


#' Internal function to translate cDNA to peptides
#'
#' @param v cDNA character vector without ambiguous bases.
#'
#' @noRd

translate_cDNA <- function(v) {
  lapply(v, function(p) {

    # protect vector length
    tryCatch(
      {
        p %>%
          Biostrings::DNAString() %>%
          Biostrings::translate(no.init.codon = TRUE) %>%
          as.character()
      },
      error = function(e) {
        return(NA)
      }
    )
  }) %>% unlist()
}
