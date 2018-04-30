## ---- garnish_summary
#' Summarize neoepitope prediction.
#'
#' Calculate neoepitope summary statistics for priority, classic, alternative, frameshift-derived, and fusion-derived neoepitopes for each sample.
#'
#' @param dt Data table. Prediction output from `garnish_affinity`.
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download an example VCF
#'    dt <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'  # extract variants
#'    antigen.garnish::garnish_variants %>%
#'
#'  # add test MHC types
#'      .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
#'                   "H-2-Kb H-2-IAd",
#'                  "HLA-A*01:47 HLA-DRB1*03:08")] %>%
#'
#'  # predict neoepitopes
#'    antigen.garnish::garnish_affinity %>%
#'
#'  # summarize predictions
#'    antigen.garnish::garnish_summary %T>%
#'    print
#'}
#' @return
#'
#'Summary table description:
#'
#' A summary data table of `dt` by `sample_id` with the following columns:
#' * **classic_neos**: classically defined neoepitopes (CDNs); defined as mutant `nmers` predicted to bind MHC with high affinity (< 50nM \eqn{IC_{50}})
#' * **classic_top_score**: sum of the top three mutant `nmer` affinity scores (see below)
#' * **alt_neos**: alternatively defined neoepitopes (ADNs); defined as mutant `nmers` predicted to bind MHC with greatly improved affinity relative to non-mutated counterparts (10x for MHC class I and 4x for MHC class II) (see below)
#' * **alt_top_score**: sum of the top three mutant `nmer` differential agretopicity indices; differential agretopicity index (DAI) is the ratio of MHC binding affinity between mutant and wt peptide (see below). If the peptide is fusion- or frameshift-derived, the binding affinity of the closest wt peptide determined by BLAST is used to calculate DAI.
#' * **priority_neos**: mutant peptides that meet both ADN and CDN criteria, and if applicable, have a fitness_score of >= 1.
#' * **fs_neos**: mutant `nmers` derived from frameshift variants predicted to bind MHC with < 500nM \eqn{IC_{50}}
#' * **fusion_neos**: mutant `nmers` derived from fusion variants predicted to bind MHC with < 500nM \eqn{IC_{50}}
#' * **nmers**: total mutant `nmers` created
#' * **predictions**: wt and mutant predictions performed
#' * **mhc_binders**: `nmers` predicted to at least minimally bind MHC (< 5000nM \eqn{IC_{50}})
#' * **fitness_scores**: Sum of the top 3 fitness_score values per sample. See `?garnish_fitness`.
#' * **garnish_score**: Sample level immune fitness. Derived by summing the exponential of `fitness_scores` for each top neoepitope across all clones in the tumor sample. See ?garnish_predicitons.
#' * **variants**: total genomic variants evaluated
#' * **transcripts**: total transcripts evaluated
#'
#'**Additional information**
#'
#'**Differential agretopicity index (DAI)** expresses the degree to which peptide binding to MHC class I or II differ due to the presence of a non-synonymous mutation. **Alternatively defined neoepitopes (ADNs)** are mutant peptides predicted to bind MHC class I or II with greatly improved affinity relative to non-mutated counterparts (*i.e.* peptides with high DAI). In mice, selection of peptides with high DAI results in a substantially improved rate of experimentally validated epitopes that mediate protection from tumor growth.
#'
#'To determine DAI, mutant peptides that at least minimally bind MHC class I or II (> 5000nM \eqn{IC_{50}}) are selected and then DAI is calculated as the fold-change in binding affinity between non-mutant and mutant peptides. ADNs are identified as mutant peptides with DAI > 10 for MHC class I and > 4 for class II.
#'
#'ADNs are generated from selective mutations in the peptide-MHC anchor position (*i.e.* the agretope) rather than mutations randomly occurring across the peptide sequence. This feature leads to two potentially important and unique immunological characteristics of ADNs. First, unlike CDNs, the TCR-facing peptide sequence in ADNs is likely the same as the corresponding non-mutant peptide. Second, the MHC binding of the corresponding non-mutant peptide may be so low that its presentation in the thymus is minimal and central tolerance may be bypassed. Functionally, there is evidence of strong immunogenicity of ADNs. Peptides, which in retrospect satisfy ADN selection criteria, are common among human tumor antigens that have been experimentally confirmed to be immunogenic. A recent extensive analysis of tumor immunity in a patient with ovarian carcinoma showed that the top five reactive mutant peptides had substantially higher mutant to non-mutant predicted MHC class I binding affinity. Moreover, a re-analysis of validated neoepitopes from non-small cell lung carcinoma or melanoma patients showed that one third of these were ADNs (resulting from an anchor position substitution that improved MHC affinity > 10-fold).
#'
#'To better model potential for oligoclonal antitumor responses directed against neoepitopes, we additionally report a **top three neoepitope score**, which is defined as the sum of the top three affinity scores \eqn{\left(\frac{1}{IC_{50}}\right)} for CDNs or sum of top three DAI for ADNs. The top three was chosen in each case because this is the minimum number that captures the potential for an oligoclonal T cell response and mirrors experimentally confirmed oligoclonality of T cell responses against human tumors. Moreover, the top three score was the least correlated to total neoepitope load (vs. top 4 through top 15) in a large scale human analysis of neoepitope across 27 disease types (R-squared = 0.0495), and therefore not purely a derivative of total neoepitope load.
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_plot}}
#' @export garnish_summary
#'
#' @references
#' Carreno BM, Magrini V, Becker-Hapak M, Kaabinejadian S, Hundal J, et al. 2015. A dendritic cell vaccine increases the breadth and diversity of melanoma neoantigen-specific T cells. Science. 348(6236):aaa3828–808
#'
#' Duan F, Duitama J, Seesi Al S, Ayres CM, Corcelli SA, et al. 2014. Genomic and bioinformatic profiling of mutational neoepitopes reveals new rules to predict anticancer immunogenicity. Journal of Experimental Medicine. 211(11):2231–48
#'
#' Fritsch EF, Rajasagi M, Ott PA, Brusic V, Hacohen N, Wu CJ. 2014. HLA-binding properties of tumor neoepitopes in humans. Cancer Immunology Research. 2(6):522–29
#'
#' Gros A, Parkhurst MR, Tran E, Pasetto A, Robbins PF, et al. 2016. Prospective identification of neoantigen-specific lymphocytes in the peripheral blood of melanoma patients. Nature Medicine. 22(4):433–38
#'
#' Jiménez-Sánchez A, Memon D, Pourpe S, Veeraraghavan H, Li Y, et al. 2017. Heterogeneous Tumor-Immune Microenvironments among Differentially Growing Metastases in an Ovarian Cancer Patient. Cell. 170(5):927–938.e20
#'
#' McGranahan N, Furness AJS, Rosenthal R, Ramskov S, Lyngaa R, et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science. 351(6280):1463–69
#'
#' Verdegaal EME, de Miranda NFCC, Visser M, Harryvan T, van Buuren MM, et al. 2016. Neoantigen landscape dynamics during human melanoma-T cell interactions. Nature. 536(7614):91–95
#'
#' @md

garnish_summary <- function(dt){

  # magrittr version check, this will not hide the error, only the NULL return on successful exit

  invisible(check_dep_versions())

# summarize over unique nmers

  dt %<>% data.table::as.data.table %>%
    data.table::copy %>%
    unique(by = c("pep_type",
                  "MHC",
                  "nmer",
                  "sample_id"))

  # set NA DAI to 0 to filter Inf and -Inf
    dt[, DAI := DAI %>% as.numeric]
    dt[is.na(DAI), DAI  := 0]


  dt <- dt[DAI != Inf &
					 DAI != -Inf &
					 Consensus_scores != Inf &
					 Consensus_scores != -Inf]

# function to sum top values of a numeric vector

sum_top_v <- function(x, value = 3){

      x %<>%
        stats::na.omit %>%
        sort %>%
        rev

      # added this in case 3 fitness scores not available, would also be useful in a sample with less than 3 nmers I guess.
      if (length(x) < value) value <- length(x)

      return(sum(x[1:value], na.rm = TRUE))
    }

  dt %>% data.table::setkey(sample_id)

dtn <- parallel::mclapply(dt[, sample_id %>% unique], function(id){

    dt <- dt[sample_id == id]

    if (!("fusion_uuid" %chin% names(dt)) %>% any) dt[, fusion_uuid := NA]
    if (!("effect_type" %chin% names(dt)) %>% any) dt[, effect_type := NA]

    if ("blast_uuid" %chin% names(dt)) dt[is.na(dai_uuid) & !is.na(blast_uuid), DAI := BLAST_A]

      return(
          data.table::data.table(
          sample_id = id,
          priority_neos_class_I = dt[class == "I" &
                                        DAI > 10 &
                                         Consensus_scores < 50] %>% nrow,
          priority_neos_class_II = dt[class == "II" &
                                        DAI > 10 &
                                        Consensus_scores < 50] %>% nrow,
          classic_neos_class_I = dt[class == "I" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 50] %>%
                                      data.table::as.data.table %>%  nrow,
          classic_neos_class_II = dt[class == "II" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 50] %>%
                                      data.table::as.data.table %>%  nrow,
          classic_top_score_class_I = dt[class == "I" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 1000, (1/Consensus_scores) %>% sum_top_v],
          classic_top_score_class_II = dt[class == "II" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 1000, (1/Consensus_scores) %>% sum_top_v],
          alt_neos_class_I = dt[class == "I" &
                                      !is.na(DAI) &
                                      Consensus_scores < 1000 &
                                      DAI > 10] %>%
                                      data.table::as.data.table %>% nrow,
          alt_neos_class_II = dt[class == "II" &
                                      !is.na(DAI) &
                                      Consensus_scores < 1000 &
                                      DAI > 10] %>%
                                      data.table::as.data.table %>%  nrow,
          alt_top_score_class_I = dt[class == "I" &
                                      !is.na(DAI) &
                                      Consensus_scores < 1000, DAI %>% sum_top_v],
          alt_top_score_class_II = dt[class == "II" &
                                      !is.na(DAI) &
                                      Consensus_scores < 1000, DAI %>% sum_top_v],
          fs_neos_class_I = dt[class == "I" &
                                      effect_type %like% "frameshift" &
                                      Consensus_scores < 1000] %>%
                                      data.table::as.data.table %>% nrow,
          fs_neos_class_II = dt[class == "II" &
                                      effect_type %like% "frameshift" &
                                      Consensus_scores < 1000]  %>%
                                      data.table::as.data.table %>% nrow,
          fusion_neos_class_I = dt[class == "I" &
                                      !is.na(fusion_uuid) &
                                      Consensus_scores < 1000] %>%
                                      data.table::as.data.table %>% nrow,
          fusion_neos_class_II = dt[class == "II" &
                                      !is.na(fusion_uuid) &
                                      Consensus_scores < 1000]  %>%
                                      data.table::as.data.table %>% nrow,
          mhc_binders_class_I = dt[class == "I" &
                                      Consensus_scores < 1000] %>% nrow,
          mhc_binders_class_II = dt[class == "II" &
                                      Consensus_scores < 1000] %>% nrow,
          predictions = dt[pep_type %like% "mut"] %>% nrow,
          nmers = dt[pep_type %like% "mut", nmer %>% unique] %>% length
          ))

    }) %>% data.table::rbindlist

# get fitness scores if available

  if ("fitness_score" %chin% names(dt)){

    for (i in dt[, sample_id %>% unique]){

      dtn[sample_id == i,
        fitness_scores_class_I :=
            dt[class == "I" & !is.na(fitness_score) & sample_id == i,
              fitness_score %>% sum_top_v]]

      dtn[sample_id == i,
        fitness_scores_class_II :=
          dt[class == "II" & !is.na(fitness_score) & sample_id == i,
            fitness_score %>% sum_top_v]]

      dtn[sample_id == i,
        priority_neos_class_I :=
          dt[class == "I" &  DAI > 10 & Consensus_scores < 50 &
          (nchar(nmer) != 9 || fitness_score >= 1) &
          sample_id == i,
            nmer_uuid %>% length]]

      dtn[sample_id == i,
        priority_neos_class_II :=
          dt[class == "II" &  DAI > 10 & Consensus_scores < 50 &
            (nchar(nmer) != 9 || fitness_score >= 1) &
              sample_id == i,
                nmer_uuid %>% length]]

      }

  }

  if ("garnish_score" %chin% names(dt))
    dtn <- merge(dtn, dt[!is.na(garnish_score), .SD %>% unique, .SDcols = c("sample_id", "garnish_score")],
      all.x = TRUE, by = "sample_id")

  if (c("ensembl_transcript_id", "var_uuid") %chin% (dt %>% names) %>% all) {

dtnv <- parallel::mclapply(dt[, sample_id %>% unique], function(id){

        dt <- dt[sample_id == id]

          return(
              data.table::data.table(
              sample_id = id,
              variants = dt[, var_uuid %>% unique] %>% length,
              transcripts = dt[ensembl_transcript_id %>% unique] %>% nrow
              ))

        }) %>% data.table::rbindlist

        dtn <- merge(dtn, dtnv, by = "sample_id", all = TRUE)
    }

    return(dtn)
}



## ---- garnish_variants
#' Process VCF variants and return a data table for epitope prediction.
#'
#' Process paired tumor-normal VCF variants annotated with [SnpEff](http://snpeff.sourceforge.net/) for neoepitope prediction using `garnish_affinity`. VCFs from matched samples can be optionally intersected to select only variants present across input files.
#'
#' @param vcfs Paths to one or more VFC files to import. See details below.
#'
#' @return A data table with one unique SnpEff variant annotation per row, including:
#' * **sample_id**: sample identifier constructed from input \code{.bam} file names
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
#'
#' @details `vcf`s can optionally contain an `AF` or `CF` *INFO* field, in which case cellular fraction or allelic fraction is considered when ranking neoepitopes. [Example vcf](http://get.rech.io/antigen.garnish_example.vcf). Single samples are required. Multi-sample `vcf`s are not supported. `vcf`s must be annotated with [SnpEff](http://snpeff.sourceforge.net/).
#'
#' Recommended workflow:
#'
#' 1. Call variants using [MuTect2](https://github.com/broadinstitute/gatk) and [Strelka2](https://github.com/Illumina/strelka), [intersecting variants](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5394620/).
#' 3. Variant filtering according to experimental specifications.
#' 4. Classification by somatic status and clonality using [PureCN](http://www.github.com/lima1/PureCN).
#' 5. Annotate variants using [SnpEff](http://snpeff.sourceforge.net/) (**required**).
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download an example VCF
#'    dt <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'  # extract variants
#'    garnish_variants %T>%
#'    str
#'}
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

garnish_variants <- function(vcfs){

  # magrittr version check, this will not hide the error, only the NULL return on successful exit

  invisible(check_dep_versions())

  message("Loading VCFs")

sdt <- lapply(vcfs %>% seq_along, function(ivf){

  # load dt
      vcf <-  vcfR::read.vcfR(vcfs[ivf], verbose = TRUE)

     if (!is.null(vcf@gt)){
      sample_id <- vcf@gt %>%
      data.table::as.data.table %>%
      names %>%
      .[-1] %>% paste(collapse = ".")
    }

  # fallback if sample names are missing

			if (sample_id == ""){
			  warning(paste0(
			  "No sample names in ", vcfs[ivf], ", using file name."))
			  sample_id <- basename(vcfs[ivf])
			  }

  # extract vcf type

      vcf_type <- vcf@meta %>%
        unlist %>%
        stringr::str_extract(stringr::regex("strelka|mutect|varscan|samtools|somaticsniper|freebayes|virmid",ignore_case = TRUE)) %>%
        stats::na.omit %>%
        unlist %>%
        data.table::first

      if (vcf_type %>% length == 0)
      	vcf_type <- "unknown"

  # return a data table of variants

    vdt <- vcf %>% get_vcf_info_dt

  # check that VCF is SnpEff-annotated

    if (
        vdt[, INFO %>% unique] %>% is.na ||
        !vdt[, INFO %>%
      stringr::str_detect(stringr::fixed("ANN=")) %>% all]
      )
      stop(paste0("\nInput file \n",
									vcfs[ivf],
									"\nis missing INFO SnpEff annotations"))

  	# parse sample level info

	    if (vcf@gt %>% length > 0)
	      vdt %<>% cbind(vcf %>% get_vcf_sample_dt)

    if (vdt %>% nrow < 1)
    	return(data.table::data.table(sample_id = sample_id))

    vdt[, sample_id := sample_id]
    vdt[, vcf_type := vcf_type]


		if (vdt %>% nrow < 1){
			warning("No variants are present in the input file.")
			return(data.table::data.table(sample_id = sample_id))
			}


		if (vdt %>% class %>%
		     .[1] == "try-error"){
		  warning("Error parsing input file INFO field.")
		  return(data.table::data.table(sample_id = sample_id))
      }

		# parse ANN column

		vdt %<>% get_vcf_snpeff_dt

		if (vdt %>% class %>%
		     .[1] == "try-error"){
		  warning("Error parsing input file SnpEff ANN annotation.")
		  return(data.table::data.table(sample_id = sample_id))
      }

    # filter SnpEff warnings

    vdt %<>% .[!ANN %like% "ERROR_.*CHROMOSOME"]
    vdt %<>% .[!ANN %likef% "WARNING_SEQUENCE_NOT_AVAILABLE"]
    vdt %<>% .[!ANN %likef% "WARNING_TRANSCRIPT_INCOMPLETE"]
    vdt %<>% .[!ANN %likef% "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS"]
    vdt %<>% .[!ANN %likef% "WARNING_TRANSCRIPT_NO_START_CODON"]

    if (vdt %>% nrow < 1){
      warning("No variants are present in the input file after filtering.")
    	return(data.table::data.table(sample_id = sample_id))
    	}

    # filter out NA
    vdt %<>% .[!ensembl_transcript_id %>% is.na &
               !cDNA_change %>% is.na]

    if (vdt %>% nrow < 1){
      warning("No variants are present in the input file after filtering.")
    	return(data.table::data.table(sample_id = sample_id))
    	}

    if ("CF" %chin% (vdt %>% names))
      vdt <- data.table::setnames("CF", "cellular_fraction")

    if ("AF" %chin% (vdt %>% names))
      vdt <- data.table::setnames("AF", "allelic_fraction")

    return(vdt)
    })

  if(class(sdt)[1] == "list")
    sdt %<>% data.table::rbindlist(use.names = TRUE,
                                   fill = TRUE)

  if (nrow(sdt) == 0){
    warning("No samples returned passing variants.")
    return(NULL)
  }

# select protein coding variants without NA
sdt %<>%
  .[protein_coding == TRUE &
  !protein_change %>% is.na &
  !effect_type %>% is.na &
  effect_type %like% "insertion|deletion|missense|frameshift"]

  if (nrow(sdt) == 0){
    warning("No samples returned protein coding variants.")
    return(NULL)
  }

  return(sdt)

}



## ---- garnish_plot
#' Graph `garnish_affinity` results.
#'
#' Plot ADN, CDN, priority, frameshift, and fusion derived `nmers` for class I and class II MHC by `sample_id`.
#'
#' @param input Output from `garnish_affinity` to graph. `input` may be a data table object, list of data tables, or file path to a rio::import-compatible file type. If a list of data tables is provided, unique plots will be generated for each data table.
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_summary}}
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download example output
#'    "antigen.garnish_example_jaffa_output.txt" %T>%
#'    utils::download.file(
#'     "http://get.rech.io/antigen.garnish_example_jaffa_output.txt", .) %>%
#'
#'  # create plot
#'    antigen.garnish::garnish_plot
#'
#'}
#'
#' @return NULL
#'
#' As a side effect: saves graph illustrating the number of neoepitopes in each sample and fitness model results to the working directory. The threshold for inclusion of fusion and frameshift-derived neoepitopes is \eqn{IC_{50}} < 1000nM.
#'
#' @export garnish_plot
#'
#' @md

garnish_plot <- function(input){

  # magrittr version check, this will not hide the error, only the NULL return on successful exit
  invisible(check_dep_versions())

  # check input
  if (class(input)[1] == "list" & class(input[[1]])[1] != "data.table")
    stop("Input must be a path to a rio::import-supported file type, a data.table, or a list of data tables (e.g. garnish_plot(list(dt1, dt2, dt3))")

  if (class(input)[1] == "character")
    input <- rio::import(input) %>% data.table::as.data.table

  if (class(input)[1] != "list" & class(input)[1] != "data.table")
    stop("Input must be a full file path to a rio::import-supported file type, a data.table object, or a list of data.tables (e.g. garnish_plot(list(dt1, dt2, dt3))")

  # set up graphing
      gplot_theme <-
        ggplot2::theme(
          line = ggplot2::element_line(colour = "#000000")) +
        ggplot2::theme(
          axis.line = ggplot2::element_line(color="black")) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = ggplot2::rel(2*1.2),
              colour = "#000000", lineheight = 0.9, face = "bold", vjust = 0)) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(colour = "#000000",
              size = ggplot2::rel(2*1.425), face = "bold")) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text(colour = "#000000",
              size = ggplot2::rel(2*1.425), face = "bold", angle = 90, vjust = 0)) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(colour = "#000000",
              size = ggplot2::rel(2*1.425), face = "bold", angle = 30, hjust = 0.9, vjust = 0.92)) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(colour = "#000000",
              size = ggplot2::rel(2*1.425), face = "bold", hjust = 0.9, vjust = 0.92, angle = 30)) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(colour = "#000000",
              size = ggplot2::rel(2*1))) +
        ggplot2::theme(
          strip.background = ggplot2::element_rect(fill = "#eeeeee")) +
        ggplot2::theme(
          strip.text = ggplot2::element_text(face = "bold")) +
        ggplot2::theme(
          legend.key = ggplot2::element_blank()) +
        ggplot2::theme(
          legend.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank()) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
        ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
        ggplot2::theme(
          plot.margin=grid::unit(c(0, 0, 0, 0 ), "cm")) +
        ggplot2::theme(
          strip.text.x = ggplot2::element_text(size = ggplot2::rel(2)))

      gplot_labs <- ggplot2::labs(x = "", y = "peptides")

      gplot_col <-   c("#ff80ab",
                       "#b388ff",
                       "#82b1ff",
                       "#ff9e80",
                       "#a7ffeb",
                       "#f4ff81",
                       "#b9f6ca",
                       "#ffe57f",
                       "#ff8a80",
                       "#ea80fc",
                       "#8c9eff",
                       "#80d8ff",
                       "#84ffff",
                       "#ccff90",
                       "#ffff8d",
                       "#ffd180")

      gplot_fn <- format(Sys.time(), "%d/%m/%y %H:%M:%OS") %>%
                    stringr::str_replace_all("[^A-Za-z0-9]", "_") %>%
                    stringr::str_replace_all("[_]+", "_")

  if(class(input)[1] != "list") input <- list(input)

  # function to shorten names for display

gplot_names <- function(gg_dt){

      if (any(gg_dt[, sample_id %>% unique %>% nchar] > 7)){

        for (i in gg_dt[nchar(sample_id) > 7, sample_id %>% unique] %>% seq_along){

          gg_dt[sample_id == gg_dt[nchar(sample_id) > 7, sample_id %>% unique][i],
                sample_id := sample_id %>%
                  substr(1, 7) %>% paste0(., "_", i)]
        }
      }

    return(gg_dt)

    }

  # function to fill in missing combinations of factors
  # to graph dt with even bars per sample_id

gplot_missing_combn <- function(dt){

      ns <- dt[, sample_id %>% unique %>% length]

type <- lapply(dt[, type %>% unique], function(x){
               replicate(ns * 2, x)
               }) %>% unlist

      gdt <- data.table(sample_id = dt[, sample_id %>% unique],
                        MHC = c(replicate(ns, "MHC Class I"),
                                replicate(ns, "MHC Class II")),
                        type = type,
                        N = 0) %>% unique
      return(gdt)

        }

lapply(input %>% seq_along, function(i){

    dt <- input[[i]]
    dt <- data.table::copy(dt) %>%
      unique(by = c("pep_type",
                    "MHC",
                    "nmer"))

      if (!(c("nmer",
              "MHC",
              "sample_id",
              "DAI",
              "frameshift",
              "Consensus_scores") %chin% names(dt)) %>% any)
        stop("'sample_id', 'nmer', 'MHC', 'frameshift', 'DAI', and 'Consensus_scores' columns are required in all inputs.")


    # create and filter data table
    dt <- dt[pep_type != "wt"] %>% unique

    if (!"dai_uuid" %chin% names(dt)) dt[, DAI := NA %>% as.numeric]

    if ("blast_uuid" %chin% names(dt))
      dt[is.na(dai_uuid) & !is.na(blast_uuid), DAI := BLAST_A]

  # add ADN, CDN, priority classification to table
    dt <- dt %>% .[Consensus_scores < 1000] %>%
      .[pep_type == "mutnfs" & Consensus_scores < 50, type := "CDN"] %>%
      .[!is.na(DAI) & DAI > 10, type := "ADN"] %>%
      .[Consensus_scores < 50 & DAI > 10, type := "priority"] %>%
      .[frameshift == TRUE & Consensus_scores < 1000, type := "frameshift"]

  # check if fusions present in input dt
    if (names(dt) %like% "fusion_uuid" %>% any)
      dt[!is.na(fusion_uuid) & fusion_uuid != "" &
         Consensus_scores < 1000,
         type := "fusion"]

    if ("fitness_score" %chin% names(dt))
      dt[Consensus_scores < 50 & DAI > 10 &
          (fitness_score >= 1 | nchar(nmer) != 9),
            type := "priority"]

    dt <- dt[!is.na(type)]

  # check if anything is left in the dt
    if (nrow(dt) < 1){
      warning(paste0("No neoeptiopes with Consensus_scores < 5000nM or that meet minimum classification criteria in input # ", i, " skipping to next input."))

      return(NULL)
    }

  # cat MHC alleles together by class for graphing
    dt[class == "I", MHC := "MHC Class I"]
    dt[class == "II", MHC := "MHC Class II"]


    gg_dt <- dt[, .N, by = c("sample_id", "MHC", "type")]

  gdt <- dt %>% gplot_missing_combn

  gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>%
        .[, N := max(N), by = c("sample_id", "type", "MHC")] %>% unique

  # shorten names for display

    gg_dt %<>% gplot_names

  # make first summary plot, neo class by sample_id and MHC

   g <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
            ggplot2::geom_col(ggplot2::aes(fill = type),
                              col = "black", position = "dodge") +
            ggplot2::facet_wrap(~MHC) +
            gplot_theme +
            gplot_labs +
            ggplot2::theme(legend.position = "bottom",
                           legend.title = ggplot2::element_blank()) +
            ggplot2::scale_fill_manual(values = gplot_col) +
            ggplot2::ggtitle(paste0("antigen.garnish summary"))

            ggplot2::ggsave(plot = g,
                  paste0("antigen.garnish_Neoepitopes_summary_",
                    gplot_fn,
                    ".pdf")
                  , height = 6, width = 9)

  if (nrow(dt[frameshift == TRUE]) > 0){

    dt_pl <- dt[frameshift == TRUE]

    dt_pl[, binding := "<1000nM"] %>%
      .[Consensus_scores < 500, binding := "<500nM"] %>%
      .[Consensus_scores < 50, binding := "<50nM"]

    gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "binding")]

    gdt <- dt_pl %>% gplot_missing_combn

    gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>%
      .[, N := max(N), by = c("sample_id", "binding", "MHC")] %>% unique

    gg_dt %<>% gplot_names

    # make frameshift summary plot, binding affinity by sample_id and MHC
      g <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
              ggplot2::geom_col(ggplot2::aes(fill = binding),
                                col = "black", position = "dodge") +
              ggplot2::facet_wrap(~MHC) +
              gplot_theme +
              gplot_labs +
              ggplot2::theme(legend.position = "bottom",
                             legend.title = ggplot2::element_blank()) +
              ggplot2::scale_fill_manual(values = gplot_col[1:3]) +
              ggplot2::ggtitle(paste0("Frameshift neoepitopes"))

      ggplot2::ggsave(plot = g, paste0("antigen.garnish_Frameshift_summary_",
                        gplot_fn,
                        ".pdf"), height = 6, width = 9)

      }

  if ("fusion" %chin% dt[, type %>% unique]%>% any){

    dt_pl <- dt[type == "fusion"]

    dt_pl[, binding := "<1000nM"] %>%
      .[Consensus_scores < 500, binding := "<500nM"] %>%
      .[Consensus_scores < 50, binding := "<50nM"]

    gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "binding")]

    gdt <- dt_pl %>% gplot_missing_combn

    gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>%
      .[, N := max(N), by = c("sample_id", "binding", "MHC")] %>% unique

    gg_dt %<>% gplot_names

    # make fusions summary plot, binding affinity by sample_id and MHC
      g <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
              ggplot2::geom_col(ggplot2::aes(fill = binding),
                                col = "black", position = "dodge") +
              ggplot2::facet_wrap(~MHC) +
              gplot_theme +
              gplot_labs +
              ggplot2::theme(legend.position = "bottom",
                             legend.title = ggplot2::element_blank()) +
              ggplot2::scale_fill_manual(values = gplot_col[1:3]) +
              ggplot2::ggtitle(paste0("Fusion neoepitopes"))

      ggplot2::ggsave(plot = g,
                      paste0("antigen.garnish_Fusions_summary_",
                      gplot_fn,
                      ".pdf"), height = 6, width = 9)
      }

    })

lapply(input %>% seq_along, function(i){

      score_dt <- input[[i]] %>% garnish_summary
      cols <- c("sample_id", names(score_dt) %include% "score(s)?_")
      score_dt <- score_dt[, .SD, .SDcols = cols]
      score_dt %<>%
      	data.table::melt(id.vars = "sample_id",
			                   measure.vars = names(score_dt) %include%
                                     "score")

      score_dt[, MHC := "MHC Class I"] %>%
        .[variable %like% "class_II", MHC := "MHC Class II"] %>%
        .[, variable := variable %>%
            stringr::str_extract("^.*(?=(_class_))")]

      if (score_dt %>% stats::na.omit %>% nrow < 1)
      	return(NULL)

      score_dt <- score_dt[!(MHC == "MHC Class II" & variable == "fitness_scores")]

      score_dt %<>% gplot_names

      if (nrow(score_dt[variable == "classic_top_score"]) != 0){

      g <- ggplot2::ggplot(score_dt[variable == "classic_top_score"],
                           ggplot2::aes(x = sample_id, y = value)) +
           ggplot2::geom_col(ggplot2::aes(fill = variable), col = "black",
                             position = "dodge") +
           ggplot2::facet_wrap(~MHC) +
           ggplot2::scale_fill_manual(values = gplot_col[1]) +
           gplot_theme +
           ggplot2::labs(x = "", y = "Score") +
           ggplot2::theme(legend.position = "none") +
           ggplot2::ggtitle(paste0("ag_classic_top_scores"))

      ggplot2::ggsave(plot = g,
        paste0("antigen.garnish_classic_scores_",
        gplot_fn,
        ".pdf"), height = 6, width = 9)

      }

      if (nrow(score_dt[variable == "alt_top_score"]) != 0){

      g <- ggplot2::ggplot(score_dt[variable == "alt_top_score"], ggplot2::aes(x = sample_id, y = value)) +
           ggplot2::geom_col(ggplot2::aes(fill = variable), col = "black",
                             position = "dodge") +
           ggplot2::facet_wrap(~MHC) +
           ggplot2::scale_fill_manual(values = gplot_col[2]) +
           gplot_theme +
           ggplot2::labs(x = "", y = "Score") +
           ggplot2::theme(legend.position = "none") +
           ggplot2::ggtitle(paste0("ag_alt_top_scores"))

      ggplot2::ggsave(plot = g,
        paste0("antigen.garnish_alt_scores_",
        gplot_fn,
        ".pdf"), height = 6, width = 9)

      }

      if (nrow(score_dt[variable == "fitness_scores"]) != 0){

      g <- ggplot2::ggplot(score_dt[variable == "fitness_scores"], ggplot2::aes(x = sample_id, y = value)) +
           ggplot2::geom_col(ggplot2::aes(fill = variable), col = "black",
                             position = "dodge") +
           ggplot2::facet_wrap(~MHC) +
           ggplot2::scale_fill_manual(values = gplot_col[3]) +
           gplot_theme +
           ggplot2::labs(x = "", y = "Score") +
           ggplot2::theme(legend.position = "none") +
           ggplot2::ggtitle(paste0("ag_fitness_scores"))

      ggplot2::ggsave(plot = g,
        paste0("antigen.garnish_fitness_summary_",
        gplot_fn,
        ".pdf"), height = 6, width = 9)
      }

  })

  lapply(input %>% seq_along, function(i){

    score_dt <- input[[i]] %>% garnish_summary

    if (!"garnish_score" %chin% names(score_dt))
    	return(NULL)

    score_dt %<>% .[, .SD %>% unique, .SDcols = c("sample_id", "garnish_score")]

    if (score_dt %>% stats::na.omit %>% nrow < 1)
    	return(NULL)

    score_dt[is.na(garnish_score), garnish_score := 0]

    score_dt %<>% gplot_names

    g <- ggplot2::ggplot(score_dt,
                         ggplot2::aes(x = sample_id, y = garnish_score)) +
         ggplot2::geom_col(fill = gplot_col[1], col = "black", position = "dodge") +
         gplot_theme +
         ggplot2::labs(x = "", y = "garnish_score") +
         ggplot2::theme(legend.position = "none") +
         ggplot2::ggtitle(paste0("antigen_garnish_score"))

    ggplot2::ggsave(plot = g,
      paste0("antigen.garnish_garnish_score_",
      gplot_fn,
      ".pdf"), height = 6, width = 9)

})

  return(NULL)

}


## ---- garnish_antigens
#' List top neoepitopes for each sample and/or by clones within each sample..
#'
#' @param dt An output data table from `garnish_affinity`.
#'
#' @return A data table with the top neoepitope per sample, and if possible per clone, in rank order of clone frequency.
#'
#' @export garnish_antigens
#' @md

garnish_antigens <- function(dt){

  if (class(dt)[1] == "character") dt <- dt %>%
  rio::import %>%
    data.table::as.data.table

  if (class(dt)[1] == "data.frame") dt %<>%
    data.table::as.data.table

  dt %<>% data.table::copy

  c <- c("clone_id", "cl_proportion")

  if (!"clone_id" %chin% names(dt)) c <- NULL

  dt[, fs := max(fitness_score, na.rm = TRUE), by = c("sample_id", c[1])]

  dt <- dt[fitness_score == fs]

  n <- names(dt)[which(names(dt) %chin% c("cDNA_change", "protein_change"))]

  if (length(n) < 1) n <- NULL

  dt <- dt[, .SD %>% unique, .SDcols = c("sample_id", "nmer", "MHC", "external_gene_name", n,
                                        "Consensus_scores", "fitness_score", "iedb_score", "min_DAI", c)]

  if (!"clone_id" %chin% names(dt)) dt <- dt %>% .[order(sample_id)]

  if ("clone_id" %chin% names(dt)) dt <- dt %>% .[order(sample_id, clone_id)]

  return(dt[Consensus_scores < 1000])

}
