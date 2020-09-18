


#' Process VCF variants and return a data table for epitope prediction.
#'
#' Process paired tumor-normal VCF variants annotated with [SnpEff](http://snpeff.sourceforge.net/) for neoantigen prediction using `garnish_affinity`. All versioned Ensembl transcript IDs (e.g. ENST00000311936.8) from any GRCh38 or GRCm38 release are supported.
#'
#' @param vcfs Paths to one or more VFC files to import. See details below.
#' @param tumor_sample_name Character, name of column in vcf of tumor sample.
#'
#' @return A data table with one unique SnpEff variant annotation per row, including:
#' * **sample_id**: sample identifier constructed from input \code{.vcf} file names
#' * **ANN**: SnpEff annotation
#' * **effect_type**: SnpEff effect type
#' * **ensembl_transcript_id**: transcript effected
#' * **cDNA_change**: cDNA change in [HGVS](http://varnomen.hgvs.org/recommendations/protein/) format
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_antigens}}
#'
#' @export garnish_variants
#' @md

garnish_variants <- function(vcfs, tumor_sample_name = "TUMOR") {
  message("Loading VCFs")

  sdt <- lapply(vcfs %>% seq_along(), function(ivf) {

    # load dt

    vcf <- vcfR::read.vcfR(vcfs[ivf],
      checkFile = FALSE,
      verbose = TRUE,
      check_keys = FALSE
    )

    sample_id <- basename(vcfs[ivf])

    # extract vcf type

    if (nrow(vcf@fix) == 0) {
      warnings("No lines in VCF.")
      return(NULL)
    }


    vcf_type <- vcf@meta %>%
      unlist() %>%
      stringr::str_extract(stringr::regex("strelka|mutect|varscan|samtools|somaticsniper|freebayes|virmid", ignore_case = TRUE)) %>%
      stats::na.omit() %>%
      unlist() %>%
      data.table::first()

    if (vcf_type %>% length() == 0) {
      vcf_type <- "unknown"
    }

    # rename generic columns to prevent downstream errors

    # return a data table of variants
    vdt <- vcf %>% get_vcf_info_dt()

    # check that VCF is SnpEff-annotated

    if (
      vdt[, INFO %>% unique()] %>% is.na() ||
        !vdt[, INFO %>%
          stringr::str_detect(stringr::fixed("ANN=")) %>% all()]
    ) {
      stop(paste0(
        "\nInput file \n",
        vcfs[ivf],
        "\nis missing INFO SnpEff annotations"
      ))
    }

    # parse sample level info

    if (vcf@gt %>% length() > 0) {
      vdt %<>% cbind(vcf %>% get_vcf_sample_dt())
    }

    if (vdt %>% nrow() < 1) {
      return(data.table::data.table(sample_id = sample_id))
    }

    vdt[, sample_id := sample_id]
    vdt[, vcf_type := vcf_type]


    if (vdt %>% nrow() < 1) {
      warning("No variants are present in the input file.")
      return(data.table::data.table(sample_id = sample_id))
    }


    if (vdt %>% class() %>%
      .[1] == "try-error") {
      warning("Error parsing input file INFO field.")
      return(data.table::data.table(sample_id = sample_id))
    }

    # parse ANN column
    vdt %<>% get_vcf_snpeff_dt

    if (vdt %>% class() %>%
      .[1] == "try-error") {
      warning("Error parsing input file SnpEff ANN annotation.")
      return(data.table::data.table(sample_id = sample_id))
    }

    # filter SnpEff warnings

    vdt %<>% .[!ANN %like% "ERROR_.*CHROMOSOME"]
    vdt %<>% .[!ANN %likef% "WARNING_SEQUENCE_NOT_AVAILABLE"]

    if (vdt %>% nrow() < 1) {
      warning("No variants are present in the input file after filtering.")
      return(data.table::data.table(sample_id = sample_id))
    }

    # filter out NA
    vdt %<>% .[!cDNA_change %>% is.na() & !is.na(ensembl_transcript_id)]

    if (vdt %>% nrow() < 1) {
      warning("No variants are present in the input file after filtering.")
      return(data.table::data.table(sample_id = sample_id))
    }

    # check tumor_sample_name is in vcf
    # AF is a GT field not an INFO field, and account for multiple alt alleles in this scenario
    if (length(tumor_sample_name) != 1 | class(tumor_sample_name)[1] != "character") {
      stop("A single tumor sample name must be provided.")
    }

    afield <- paste(tumor_sample_name, "_AF", sep = "")

    if (!afield %chin% (vdt %>% names())) {
      pns <- names(vdt) %include% "_AF" %>% stringr::str_replace("_AF", "")
    }

    if (afield %chin% (vdt %>% names())) {
      vdt %<>% tidyr::separate_rows(eval(afield), sep = ",")
    }

    vdt %<>% tidyr::separate_rows("ALT", sep = ",")

    vdt %<>% data.table::as.data.table(.)

    # now keep only rows that match the previously split ANN field
    vdt %<>% .[ALT == stringr::str_extract(ANN, pattern = "(<?=ANN=)?[AGCT]+(?=\\|)")]

    return(vdt)
  })

  if (class(sdt)[1] == "list") {
    sdt %<>% data.table::rbindlist(
      use.names = TRUE,
      fill = TRUE
    )
  }

  if (nrow(sdt) == 0 || ncol(sdt) == 1) {
    warning("No samples returned passing variants.")
    return(NULL)
  }

  # select protein coding variants without NA
  sdt %<>%
    .[
      !protein_change %>% is.na() &
        !effect_type %>% is.na() &
        effect_type %like% "insertion|deletion|missense|frameshift"
    ]


  if (nrow(sdt) == 0) {
    warning("No samples returned protein coding variants.")
    return(NULL)
  }

  message("\nDone loading variants.")
  return(sdt)
}



#' List top peptides using TESLA criteria for recognition features of immunogenic peptides.
#'
#' The \href{https://www.parkerici.org/research-project/tumor-neoantigen-selection-alliance-tesla/}{TESLA consortium} identified recognition features of immunogenic peptides. This function filters peptides meeting any of these criteria.
#'
#' @param dt An output data table from `garnish_affinity`, either a data table object or path to a file.
#' @param affinity_threshold Numeric. Neoantigen affinity threshold, nanomolar (nM) scale.
#' @param differential_agretopcity_threshold Numeric. Neoantigen differential agretopcity threshold. Differential agretopicty is the proteome-wide ratio of MHC binding afinity between mutant and closest normal peptide, with higher values indicating greater relative binding of the mutant peptide.
#' @param foreignness_threshold Numeric. Neoantigen foreignness threshold. Value of 0 to 1 indicating the TCR recognition probability, calculated by summing alignments in IEDB immunogenic peptides, with 1 indicating greater homology to immunogenic peptides.
#' @param dissimilarity_threshold Numeric. Neoantigen dissimilarity threshold. Value of 0 to 1 indicating alignment to the self-proteome, calculated in an analogous manner to neoanigen foreignness, with 1 indicating greater dissimilarity.
#'
#' @return A data table with ranked and annotated peptides.
#'
#' @seealso \code{\link{garnish_variants}}
#' @seealso \code{\link{garnish_affinity}}
#'
#' @references
#' Richman LP, Vonderheide RH, and Rech AJ. Neoantigen dissimilarity to the self-proteome predicts immunogenicity and response to immune checkpoint blockade. Cell Systems. 2019.

#' Duan, F., Duitama, J., Seesi, S.A., Ayres, C.M., Corcelli, S.A., Pawashe, A.P., Blanchard, T., McMahon, D., Sidney, J., Sette, A., et al. Genomic and bioinformatic profiling of mutational neoepitopes reveals new rules to predict anticancer immunogenicity. J Exp Med. 2014.
#'
#' Luksza, M, Riaz, N, Makarov, V, Balachandran VP, et al. A neoepitope fitness model predicts tumour response to checkpoint blockade immunotherapy. Nature. 2017.
#' Rech AJ, Balli D, Mantero A, Ishwaran H, Nathanson KL, Stanger BZ, Vonderheide RH. Tumor immunity and survival as a function of alternative neopeptides in human cancer. Clinical Cancer Research, 2018.
#'
#' Wells DK, van Buuren MM, Dang KK, Hubbard-Lucey VM, Sheehan KCF, Campbell KM, Lamb A, Ward JP, Sidney J, Blazquez AB, Rech AJ, Zaretsky JM, Comin-Anduix B, Ng AHC, Chour W, Yu TV, Rizvi1 H, Chen JM, Manning P, Steiner GM, Doan XC, The TESLA Consortium, Merghoub T, Guinney J, Kolom A, Selinsky C, Ribas A, Hellmann MD, Hacohen N, Sette A, Heath JR, Bhardwaj N, Ramsdell F, Schreiber RD, Schumacher TN, Kvistborg P, Defranoux N. Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction. Cell. 2020.
#' @export garnish_antigens
#' @md

garnish_antigens <- function(dt,
                             affinity_threshold = 34,
                             differential_agretopcity_threshold = 10,
                             dissimilarity_threshold = 0,
                             foreignness_threshold = 10e-16) {
  if (class(dt)[1] == "character") {
    dt <- dt %>%
      data.table::fread()
  }

  if (class(dt)[1] == "data.frame") {
    dt %<>%
      data.table::as.data.table()
  }

  dt %<>% data.table::copy()

  if (!"Ensemble_score" %chin% names(dt)) {
    message("Missing Ensemble_score column. Input to garnish_antigens must be garnish_affinity output.")
    return(dt)
  }

  dt <- dt[Ensemble_score < 34 & pep_type != "wt"]

  if (nrow(dt) == 0) {
    ("No qualifying neoantigens present in data.")
  }

  if (!"min_DAI" %chin% names(dt) & "DAI" %chin% names(dt)) {
    dt[, min_DAI := DAI]
  }

  ol <- c("dissimilarity", "foreignness_score", "min_DAI", "Ensemble_score")

  ml <- ol[which(!ol %chin% names(dt))]

  ol <- ol[which(ol %chin% names(dt))]

  message(
    paste("Ranking peptides based on the follow metrics from available input: ",
      paste(ol, collapse = ", "), ".",
      sep = ""
    )
  )

  if (length(ml) != 0) dt[, as.character(ml) := as.numeric(NA)]

  n <- names(dt)[which(names(dt) %chin%
    c("cDNA_change", "protein_change", "ensembl_transcript_id"))]

  if (length(n) < 1) n <- NULL

  dt <- dt[, .SD %>% unique(), .SDcols = c(
    "sample_id", "nmer", "MHC", n,
    "Ensemble_score", "dissimilarity", "foreignness_score", "min_DAI"
  )]

  annotate_antigens <- function(ie, md, diss) {
    iedb_l <- ifelse(ie > 10e-16, "foreignness", as.character(NA))

    dai_l <- ifelse(md > 10, "agretopicity", as.character(NA))

    diss_l <- ifelse(diss > 0, "dissimilarity", as.character(NA))

    anno <- paste("binding affinity", iedb_l, dai_l, diss_l, sep = "; ")

    anno %<>% stringr::str_replace_all("NA;|NA$", "")

    anno %<>% stringr::str_replace_all("[\\ ]+", " ")

    anno %<>% stringr::str_replace_all(";[\\ ]+?$", "")

    return(anno)
  }

  dt[, Recognition_Features := annotate_antigens(
    foreignness_score,
    min_DAI,
    dissimilarity
  )]

  return(dt)
}
