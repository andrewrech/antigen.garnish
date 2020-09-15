


#' Process VCF variants and return a data table for epitope prediction.
#'
#' Process paired tumor-normal VCF variants annotated with [SnpEff](http://snpeff.sourceforge.net/) for neoantigen prediction using `garnish_affinity`. All versioned Ensembl transcript IDs (e.g. ENST00000311936.8) from any GRCh38 or GRCm38 release are supported.
#'
#' @param vcfs Paths to one or more VFC files to import. See details below.
#' @param tumor_sample_name Character, name of column in vcf of tumor sample, used to determine mutant allelic fraction of neoantigens.
#'
#' @return A data table with one unique SnpEff variant annotation per row, including:
#' * **sample_id**: sample identifier constructed from input \code{.vcf} file names
#' * **se**: SnpEff annotation
#' * **effect_type**: SnpEff effect type
#' * **ensembl_transcript_id**: transcript effected
#' * **ensembl_gene_id**: gene effected
#' * **protein_change**: protein change in [HGVS](http://varnomen.hgvs.org/recommendations/DNA/) format
#' * **cDNA_change**: cDNA change in [HGVS](http://varnomen.hgvs.org/recommendations/protein/) format
#' * **protein_coding**: is the variant protein coding?
#'
#' if CF or AF fields in provided in input VCFs, either:
#' * **cellular_fraction**: cell fraction taken from input, such as from clonality estimates from [PureCN](http://www.github.com/lima1/PureCN)
#' * **allelic_fraction**: allelic fraction taken from input
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_antigens}}
#'
#' @details `vcf`s must be annotated with [SnpEff](http://snpeff.sourceforge.net/). `vcf`s can optionally contain an `AF` or `CF` *INFO* field, in which case cellular fraction or allelic fraction is considered when ranking neoantigens. See [example vcf](http://get.rech.io/antigen.garnish_example.vcf). Single samples are required. Multi-sample `vcf`s are not supported.
#'
#' Recommended workflow:
#'
#' 1. Call variants using [MuTect2](https://github.com/broadinstitute/gatk) and [Strelka2](https://github.com/Illumina/strelka), [intersecting variants](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5394620/).
#' 3. Variant filtering according to experimental specifications.
#' 4. Classification by somatic status and clonality using [PureCN](http://www.github.com/lima1/PureCN).
#' 5. Annotate variants using [SnpEff](http://snpeff.sourceforge.net/) (**required**).
#'
#' @examples
#' \dontrun{
#'
#' # load an example VCF
#' dir <- system.file(package = "antigen.garnish") %>%
#'   file.path(., "extdata/testdata")
#'
#' dt <- "antigen.garnish_example.vcf" %>%
#'   file.path(dir, .) %>%
#'
#'   # extract variants
#'   garnish_variants() %T>%
#'   str
#' }
#'
#' @references
#' Krøigård AB, Thomassen M, Lænkholm A-V, Kruse TA, Larsen MJ. 2016. Evaluation of Nine Somatic Variant Callers for Detection of Somatic Mutations in Exome and Targeted Deep Sequencing Data. PLoS ONE. 11(3):e0151664
#'
#' Cingolani P, Platts A, Wang LL, Coon M, Nguyen T, et al. 2012. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly (Austin). 6(2):80–92
#'
#' Riester M, Singh AP, Brannon AR, Yu K, Campbell CD, et al. 2016. PureCN: copy number calling and SNV classification using targeted short read sequencing. Source Code Biol Med. 11(1):13
#'
#' Callari M, Sammut SJ, De Mattos-Arruda L, Bruna A, Rueda OM, Chin SF, and Caldas C. 2017. Intersect-then-combine approach: improving the performance of somatic variant calling in whole exome sequencing data using multiple aligners and callers. Genome Medicine. 9:35.
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

    # return a data table of variants

    vdt <- vcf %>% get_vcf_info_dt()

    # rename generic columns to prevent downstream errors

    if (names(vdt) %like% "^V[0-9]+$" %>% any()) {
      vdt %>% data.table::setnames(
        names(vdt) %include% "^V[0-9]+$",
        paste(names(vdt) %include% "^V[0-9]+$", ".x", sep = "")
      )
    }

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
    vdt %<>% .[!ANN %likef% "WARNING_TRANSCRIPT_INCOMPLETE"]
    vdt %<>% .[!ANN %likef% "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS"]
    vdt %<>% .[!ANN %likef% "WARNING_TRANSCRIPT_NO_START_CODON"]

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

    if ("CF" %chin% (vdt %>% names())) {
      vdt %>% data.table::setnames("CF", "cellular_fraction")
    }

    # check tumor_sample_name is in vcf
    # AF is a GT field not an INFO field, and account for multiple alt alleles in this scenario
    if (length(tumor_sample_name) != 1 | class(tumor_sample_name)[1] != "character") {
      stop("A single tumor sample name must be provided.")
    }

    afield <- paste(tumor_sample_name, "_AF", sep = "")

    if (!afield %chin% (vdt %>% names())) {
      pns <- names(vdt) %include% "_AF" %>% stringr::str_replace("_AF", "")


      # if sample names exist, report mismatches
      if (!(pns == "") %>% all()) {
        message("Unable to match tumor sample name, allelic_fraction will not be considered in downstream analysis.")
        message(paste("Possible sample names from the VCF are:", paste(pns, collapse = "\n"), sep = "\n"))
        message("To include allelic_fraction in downstream analysis, rerun garnish_variants with the correct tumor_sample_name argument.")
      }
    }

    if (afield %chin% (vdt %>% names())) {
      vdt %>% data.table::setnames(afield, "allelic_fraction")

      vdt %<>% tidyr::separate_rows(c("ALT", "allelic_fraction"), sep = ",")

      vdt %<>% data.table::as.data.table(.)

      # now keep only rows that match the previously split ANN field

      vdt <- vdt[ALT == stringr::str_extract(ANN, pattern = "^[AGCT]+(?=\\|)")]
    }

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
    .[protein_coding == TRUE &
      !protein_change %>% is.na() &
      !effect_type %>% is.na() &
      effect_type %like% "insertion|deletion|missense|frameshift"]

  if (nrow(sdt) == 0) {
    warning("No samples returned protein coding variants.")
    return(NULL)
  }

  if ("cellular_fraction" %chin% names(sdt)) {
    sdt[, cellular_fraction := cellular_fraction %>% as.numeric()]
  }

  if ("allelic_fraction" %chin% names(sdt)) {
    sdt[, allelic_fraction := allelic_fraction %>% as.numeric()]
  }


  message("\nDone loading variants.")
  return(sdt)
}





#' List top neoantigens for each sample and/or by clones within each sample using TESLA criteria for recognition features of immunogenic peptides.
#'
#' The \href{https://www.parkerici.org/research-project/tumor-neoantigen-selection-alliance-tesla/}{TESLA consortium} identified recognition features of immunogenic peptides. This function applies each of those criteria and returns any neoantigen that meets MHC binding affinity (< 34nM) and annotates in the Recognition_Features column with: agretopicity (DAI > 10), foreignness (iedb_score > 10e-16), and dissimilarity > 0.
#'
#' @param dt An output data table from `garnish_affinity`, either a data table object or path to a file.
#'
#' @return A data table with all neoantigens from the input table that meet at least one of the recognition features of immunogenic peptides from the \href{https://www.parkerici.org/research-project/tumor-neoantigen-selection-alliance-tesla/}{TESLA consortium}.
#'
#' @seealso \code{\link{garnish_variants}}
#' @seealso \code{\link{garnish_affinity}}
#'
#' @export garnish_antigens
#' @md

garnish_antigens <- function(dt) {
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
    stop("Missing Ensemble_score column.  Input to garnish_antigens must be garnish_affinity output.")
  }

  dt <- dt[Ensemble_score < 34 & pep_type != "wt"]

  if (nrow(dt) == 0) stop("No qualifying neoantigens present in data.")

  if (!"min_DAI" %chin% names(dt) & "DAI" %chin% names(dt)) {
    dt[, min_DAI := DAI]
  }

  ol <- c("dissimilarity", "iedb_score", "min_DAI", "Ensemble_score")

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
    c("cDNA_change", "protein_change", "ensembl_transcript_id", "clone_id"))]

  if (length(n) < 1) n <- NULL

  dt <- dt[, .SD %>% unique(), .SDcols = c(
    "sample_id", "nmer", "MHC", n,
    "Ensemble_score", "dissimilarity", "iedb_score", "min_DAI"
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
    iedb_score,
    min_DAI,
    dissimilarity
  )]

  if (!"clone_id" %chin% names(dt)) dt <- dt %>% .[order(sample_id)]

  if ("clone_id" %chin% names(dt)) dt <- dt %>% .[order(sample_id, clone_id)]

  return(dt)
}
