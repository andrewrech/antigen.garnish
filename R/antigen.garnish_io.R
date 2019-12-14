
#' Summarize neoantigen prediction.
#'
#' Calculate neoantigen summary statistics for priority, classic, alternative, frameshift-derived, and fusion-derived neoantigens for each sample.
#'
#' @param dt Data table. Prediction output from `garnish_affinity`.
#'
#' @examples
#'\dontrun{
#'
#'# load an example VCF
#'	dir <- system.file(package = "antigen.garnish") %>%
#'		file.path(., "extdata/testdata")
#'
#'	dt <- "antigen.garnish_example.vcf" %>%
#'	file.path(dir, .) %>%
#'
#'# extract variants
#'    garnish_variants %>%
#'
#'# add space separated MHC types
#' # see list_mhc() for nomenclature of supported alleles
#'
#'     .[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")] %>%
#'
#'# predict neoantigens
#'    garnish_affinity
#'
#'# summarize predictions
#'   dt %>%
#'     garnish_summary %T>%
#'       print
#'
#'# generate summary graphs
#'   dt %>% garnish_plot
#'}
#' @return
#'
#'Summary table description:
#'
#' A summary data table of `dt` by `sample_id` with the following columns:
#' * **classic_neos**: classically defined neoantigens (CDNs); defined as mutant `nmers` predicted to bind MHC with high affinity (< 50nM \eqn{IC_{50}})
#' * **classic_top_score**: sum of the top three mutant `nmer` affinity scores (see below)
#' * **alt_neos**: alternatively defined neoantigens (ADNs); defined as mutant `nmers` predicted to bind MHC with greatly improved affinity relative to non-mutated counterparts (10x for MHC class I and 4x for MHC class II) (see below)
#' * **alt_top_score**: sum of the top three mutant `nmer` differential agretopicity indices; differential agretopicity index (DAI) is the ratio of MHC binding affinity between mutant and wt peptide (see below). If the peptide is fusion- or frameshift-derived, the binding affinity of the closest wt peptide determined by BLAST is used to calculate DAI.
#' * **fs_neos**: mutant `nmers` derived from frameshift variants predicted to bind MHC with < 500nM \eqn{IC_{50}}
#' * **fusion_neos**: mutant `nmers` derived from fusion variants predicted to bind MHC with < 500nM \eqn{IC_{50}}
#' * **IEDB_high**: defined as mutant `nmers` predicted to bind MHC with affinity (< 500nM \eqn{IC_{50}}) and IEDB alignment score > 0.9.
#' * **high_dissimilarity**: defined as mutant `nmers` predicted to bind MHC with affinity (< 500nM \eqn{IC_{50}}) and dissimilarity > 0.75.
#' * **nmers**: total mutant `nmers` created
#' * **predictions**: wt and mutant predictions performed
#' * **mhc_binders**: `nmers` predicted to at least minimally bind MHC (< 500nM \eqn{IC_{50}})
#' * **variants**: total genomic variants evaluated
#' * **transcripts**: total transcripts evaluated
#'
#'**Additional information**
#'
#'**Differential agretopicity index (DAI)** expresses the degree to which peptide binding to MHC class I or II differ due to the presence of a non-synonymous mutation. **Alternatively defined neoantigens (ADNs)** are mutant peptides predicted to bind MHC class I or II with greatly improved affinity relative to non-mutated counterparts (*i.e.* peptides with high DAI). In mice, selection of peptides with high DAI results in a substantially improved rate of experimentally validated epitopes that mediate protection from tumor growth.
#'
#'To determine DAI, mutant peptides that at least minimally bind MHC class I or II (> 500nM \eqn{IC_{50}}) are selected and then DAI is calculated as the fold-change in binding affinity between non-mutant and mutant peptides. ADNs are identified as mutant peptides with DAI > 10 for MHC class I and > 4 for class II.
#'
#'ADNs are generated from selective mutations in the peptide-MHC anchor position (*i.e.* the agretope) rather than mutations randomly occurring across the peptide sequence. This feature leads to two potentially important and unique immunological characteristics of ADNs. First, unlike CDNs, the TCR-facing peptide sequence in ADNs is likely the same as the corresponding non-mutant peptide. Second, the MHC binding of the corresponding non-mutant peptide may be so low that its presentation in the thymus is minimal and central tolerance may be bypassed. Functionally, there is evidence of strong immunogenicity of ADNs. Peptides, which in retrospect satisfy ADN selection criteria, are common among human tumor antigens that have been experimentally confirmed to be immunogenic. A recent extensive analysis of tumor immunity in a patient with ovarian carcinoma showed that the top five reactive mutant peptides had substantially higher mutant to non-mutant predicted MHC class I binding affinity. Moreover, a re-analysis of validated neoantigens from non-small cell lung carcinoma or melanoma patients showed that one third of these were ADNs (resulting from an anchor position substitution that improved MHC affinity > 10-fold).
#'
#'To better model potential for oligoclonal antitumor responses directed against neoantigens, we additionally report a **top three neoantigen score**, which is defined as the sum of the top three affinity scores \eqn{\left(\frac{1}{IC_{50}}\right)} for CDNs or sum of top three DAI for ADNs. The top three was chosen in each case because this is the minimum number that captures the potential for an oligoclonal T cell response and mirrors experimentally confirmed oligoclonality of T cell responses against human tumors. Moreover, the top three score was the least correlated to total neoantigen load (vs. top 4 through top 15) in a large scale human analysis of neoantigen across 27 disease types (R-squared = 0.0495), and therefore not purely a derivative of total neoantigen load.
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_plot}}
#' @seealso \code{\link{garnish_antigens}}
#'
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
           Ensemble_score != Inf &
           Ensemble_score != -Inf]

# function to sum top values of a numeric vector

sum_top_v <- function(x, value = 3){

      x %<>%
        stats::na.omit %>%
        sort %>%
        rev

      # edge case for sample with less than 3 nmers.
      if (length(x) < value) value <- length(x)

      return(sum(x[1:value], na.rm = TRUE))
    }

  dt %>% data.table::setkey(sample_id)

dtn <- lapply(dt[, sample_id %>% unique], function(id){

    dt <- dt[sample_id == id]

    if (!("fusion_uuid" %chin% names(dt)) %>% any) dt[, fusion_uuid := NA]
    if (!("effect_type" %chin% names(dt)) %>% any) dt[, effect_type := NA]

    if (!"min_DAI" %chin% names(dt)){

      dt[!is.na(DAI), min_DAI := DAI]

      message("Proteome-wide minimum differential agretopicity is not available from these data, defaulting to local to predict ADNs.")

    }

      return(
          data.table::data.table(
          sample_id = id,
          CDNs_class_I = dt[class == "I" &
                                      pep_type != "wt" &
                                      Ensemble_score < 50,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          CDNs_class_II = dt[class == "II" &
                                      pep_type != "wt" &
                                      Ensemble_score < 50,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          CDNs_top_score_class_I = dt[class == "I" &
                                      pep_type != "wt" &
                                      Ensemble_score < 50, (1/Ensemble_score) %>% sum_top_v],
          CDNs_top_score_class_II = dt[class == "II" &
                                      pep_type != "wt" &
                                      Ensemble_score < 50, (1/Ensemble_score) %>% sum_top_v],
          ADNs_class_I = dt[class == "I" &
                                      pep_type != "wt" &
                                      Ensemble_score < 500 &
                                      min_DAI > 10,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          ADNs_class_II = dt[class == "II" &
                                      pep_type != "wt" &
                                      Ensemble_score < 500 &
                                      min_DAI > 10,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          ADNs_top_score_class_I = dt[class == "I" &
                                      pep_type != "wt" &
                                      !is.na(min_DAI) &
                                      Ensemble_score < 500, min_DAI %>% sum_top_v],
         ADNs_top_score_class_II = dt[class == "II" &
                                      pep_type != "wt" &
                                      !is.na(min_DAI) &
                                      Ensemble_score < 500, min_DAI %>% sum_top_v],
          fs_neos_class_I = dt[class == "I" &
                                      effect_type %like% "frameshift" &
                                      Ensemble_score < 500,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          fs_neos_class_II = dt[class == "II" &
                                      effect_type %like% "frameshift" &
                                      Ensemble_score < 500,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          fusion_neos_class_I = dt[class == "I" &
                                      !is.na(fusion_uuid) &
                                      Ensemble_score < 500,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          fusion_neos_class_II = dt[class == "II" &
                                      !is.na(fusion_uuid) &
                                      Ensemble_score < 500,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          MHC_binders_class_I = dt[class == "I" & pep_type != "wt" &
                                      Ensemble_score < 500,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          MHC_binders_class_II = dt[class == "II" & pep_type != "wt" &
                                      Ensemble_score < 500,
                                      paste(nmer, MHC, sep = "_") %>% unique %>% length],
          predictions = dt[pep_type %like% "mut", paste(nmer, MHC, sep = "_") %>% unique %>% length],
          nmers = dt[pep_type %like% "mut", nmer %>% unique %>% length]
          ))

    }) %>% data.table::rbindlist

    if ("iedb_score" %chin% names(dt)){

      ie_I <- dt[class == "I" &
                 pep_type != "wt" &
                 Ensemble_score < 500 &
                 iedb_score > 0.9,
                 paste(nmer, MHC, sep = "_") %>% unique %>% length,
                 by = "sample_id"] %>%
                 data.table::setnames("V1", "IEDB_high_class_I")

      ie_II <- dt[class == "II" &
                  pep_type != "wt" &
                  Ensemble_score < 500 &
                  iedb_score > 0.9,
                  paste(nmer, MHC, sep = "_") %>% unique %>% length,
                  by = "sample_id"] %>%
                  data.table::setnames("V1", "IEDB_high_class_II")

      ie <- merge(ie_I, ie_II, by = "sample_id", all = TRUE)

      dtn <- merge(dtn, ie, by = "sample_id", all.x = TRUE)

    }

    if ("dissimilarity" %chin% names(dt)){

      st_I <- dt[class == "I" &
                 pep_type != "wt" &
                 Ensemble_score < 500 &
                 dissimilarity > 0.75,
                 paste(nmer, MHC, sep = "_") %>% unique %>% length,
                 by = "sample_id"] %>%
                 data.table::setnames("V1", "high_dissimilarity_class_I")

      st_II <- dt[class == "II" &
                  pep_type != "wt" &
                  Ensemble_score < 500 &
                  dissimilarity > 0.75,
                  paste(nmer, MHC, sep = "_") %>% unique %>% length,
                  by = "sample_id"] %>%
                  data.table::setnames("V1", "high_dissimilarity_class_II")

      st <- merge(st_I, st_II, by = "sample_id", all = TRUE)

      dtn <- merge(dtn, st, by = "sample_id", all.x = TRUE)

    }

  if (c("ensembl_transcript_id", "var_uuid") %chin% (dt %>% names) %>% all) {

dtnv <- lapply(dt[, sample_id %>% unique], function(id){

        dt <- dt[sample_id == id]

          return(
              data.table::data.table(
              sample_id = id,
              variants = dt[, var_uuid %>% unique %>% length],
              transcripts = dt[, ensembl_transcript_id %>% unique %>% length]
              ))

        }) %>% data.table::rbindlist

        dtn <- merge(dtn, dtnv, by = "sample_id", all = TRUE)
    }

    # convert appropriate NA values to 0
    NA_to_0 <- function(v){

      if (!class(v)[1] %chin% c("numeric", "integer"))
       return(v)

      v[which(is.na(v))] <- 0

      return(v)

    }

    dtn[, names(dtn) := lapply(.SD, NA_to_0), .SDcols = names(dtn)]

    return(dtn)
}



#' Return minimal neoantigen prediction information for all peptides.
#'
#' This function reduces \code{\link{garnish_affinity}} output to a more manageable table, dropping columns with less critical information and removing wild-type peptide rows.
#' For sample-level summary statistics, see \code{\link{garnish_summary}}, and for highest priority epitopes per sample, see \code{\link{garnish_antigens}}.
#'
#' @param dt Input data.table from `garnish_affinity` output.
#' @param slimmer Logical. Default is TRUE, set to false to return percentile/ranks of predictions from prediction tools.
#' @return A slimmed down data table from `garnish_affinity` output with minimal information for each peptide prediction. Contains at most the following columns:
#' * **sample_id**
#' * **external_gene_name**
#' * **protein_change**
#' * **nmer**
#' * **MHC**
#' * **class**: MHC class I or II
#' * **mhcflurry_prediction**
#' * **mhcflurry_prediction_percentile**
#' * **mhcnuggets_pred_lstm**
#' * **mhcnuggets_pred_gru**
#' * **affinity(nM)_netMHC**
#' * **affinity(nM)_netMHCpan**
#' * **\%Rank_netMHC**
#' * **\%Rank_netMHCpan**
#' * **Ensemble_score**
#' * **DAI**: Differential agretopicity of variant to corresponding wild-type, see `garnish_summary`.
#' * **min_DAI**: The most conservative DAI value based on a global alignment to the wild-type proteome.
#' * **iedb_score**
#' * **dissimilarity**
#' * **cl_proportion**: clonal cell fraction composed by the variant
#' * **antigen.garnish_class**: epitope class, one of "Classic", "Alternative", "IEDB high", "high dissimilarity", see `?garnish_summary` for thresholds.
#'
#' See `garnish_affinity` for more column descriptions.
#'
#' @examples
#'\dontrun{
#'
#'# load an example VCF
#'	dir <- system.file(package = "antigen.garnish") %>%
#'		file.path(., "extdata/testdata")
#'
#'	dt <- "antigen.garnish_example.vcf" %>%
#'	file.path(dir, .) %>%
#'
#'# extract variants
#'    garnish_variants %>%
#'
#'# add space separated MHC types
#' # see list_mhc() for nomenclature of supported alleles
#'
#'     .[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")] %>%
#'
#'# predict neoantigens
#'    garnish_affinity
#'
#'# summarize predictions
#'   dt %>%
#'     garnish_summary %T>%
#'       print
#'
#'# generate summary graphs
#'   dt %>% garnish_plot
#'}
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_summary}}
#' @seealso \code{\link{garnish_antigens}}
#'
#' @export garnish_slim
#' @md

garnish_slim <- function(dt, slimmer = TRUE){

  # magrittr version check, this will not hide the error, only the NULL return on succesful exit
  invisible(check_dep_versions())

  if (!"data.table" %chin% class(dt)[1])
    stop("Input must be a data table.")

  min_cols <- c("sample_id", "nmer", "MHC", "Ensemble_score", "class")

  if (!all(min_cols %chin% names(dt)))
    stop(paste("Table must at a minimum contain the columns:", paste(min_cols, collapse = " ")))

  # don't damage input table and operate on mutant rows only for speed
  dt <- dt %>% data.table::copy %>% .[pep_type != "wt"]

  addl_cols <- c("external_gene_name", "protein_change",
                  "DAI", "min_DAI", "cl_proportion", "iedb_score", "dissimilarity")

  addl_cols <- addl_cols[which(addl_cols %chin% names(dt))]

  # selected columns for slim MHC Class I predictionas
  class_I_prediction_cols <- c("mhcflurry_prediction", "mhcflurry_prediction_percentile",
                              "mhcnuggets_pred_lstm", "mhcnuggets_pred_gru",
                              "affinity(nM)_netMHC", "%Rank_netMHC",
                              "affinity(nM)_netMHCpan", "%Rank_netMHCpan")

  if (slimmer)
    class_I_prediction_cols %<>% .[which(!class_I_prediction_cols %like% "(%Ran)|(percentile)")]

  # selected columns for slim MHC Class II predictionas
  netmhcIIpan_cols <- c("affinity(nM)_netMHCIIpan", "%Rank_netMHCIIpan")
  netmhcII_cols <- c("affinity(nM)_netMHCII", "%Random_netMHCII")

  # determine which MHC Class II tools were used (depends on input MHC Class II alleles)
  if (("netMHCIIpan" %chin% (dt %>% names)) & ("netMHCII" %chin% (dt %>% names))){
    class_II_prediction_cols <- c(netmhcIIpan_cols, netmhcII_cols)
  } else if (("netMHCIIpan" %chin% (dt %>% names)) & !("netMHCII" %chin% (dt %>% names))){
    class_II_prediction_cols <- netmhcIIpan_cols
  } else if (!("netMHCIIpan" %chin% (dt %>% names)) & ("netMHCII" %chin% (dt %>% names))){
    class_II_prediction_cols <- netmhcII_cols
  }

  if (slimmer) class_II_prediction_cols %<>% .[which(!class_II_prediction_cols %like% "(%Ran)|(percentile)")]

  # determine if Class I, Class II or both predictions are present
  if ("class" %chin%  (dt %>% names)) {
    if ((dt[class == "I"] %>% nrow) > 0 & (dt[class == "II"] %>% nrow) > 0) {
      prediction_columns <- c(class_I_prediction_cols, class_II_prediction_cols)
    } else if ((dt[class == "I"] %>% nrow) > 0 & (dt[class == "II"] %>% nrow) == 0) {
      prediction_columns <- class_I_prediction_cols
    } else if ((dt[class == "I"] %>% nrow) == 0 & (dt[class == "II"] %>% nrow) > 0) {
      prediction_columns <- class_II_prediction_cols
    } else {
      stop("Input dt is missing MHC Class I and Class II predictions.")
    }
  }

  # classify as many neos as possible

  if (c("DAI", "min_DAI") %chin% names(dt) %>% any){

    if ("min_DAI" %chin% names(dt))
      dt[min_DAI > 10 & Ensemble_score < 500, ag_class1 := "Alternative"]
    if ("DAI" %chin% names(dt) & !"min_DAI" %chin% names(dt)){
      message("Using local differential agretopicty, proteome-wide was not provided in input table.")
      dt[DAI > 10 & Ensemble_score < 500, ag_class1 := "Alternative"]
    }

  }

  dt[Ensemble_score < 50, ag_class2 := "Classic"]

  if ("iedb_score" %chin% names(dt))
    dt[iedb_score > 0.9 & Ensemble_score < 500, ag_class3 := "IEDB high"]
  if ("dissimilarity" %chin% names(dt))
      dt[dissimilarity > 0.75 & Ensemble_score < 500, ag_class4 := "high dissimilarity"]

  n <- names(dt)[names(dt) %like% "ag_class[1-4]$"]

  dt[, antigen.garnish_class := paste(.SD %>% unlist %>% na.omit, collapse = ", "), .SDcols = n, by = 1:nrow(dt)]

  dt[antigen.garnish_class == "", antigen.garnish_class := "Unclassified"]

  cols_to_keep <- c(min_cols,
                    prediction_columns,
                    addl_cols,
                    "antigen.garnish_class")

  if (!all((cols_to_keep) %chin% names(dt))) {

    missing_name <- setdiff(cols_to_keep, names(dt))

    warning("Missing ", paste(missing_name, sep = ", "), " column(s) from input neoantigen prediction table. Returning without this column.")

    cols_to_keep <- cols_to_keep[which(!cols_to_keep %chin% missing_name)]

  }

  dt_slim <- dt[, .SD, .SDcols = cols_to_keep] %>% .[order(sample_id, Ensemble_score)]

  return(dt_slim)

}




#' Process VCF variants and return a data table for epitope prediction.
#'
#' Process paired tumor-normal VCF variants annotated with [SnpEff](http://snpeff.sourceforge.net/) for neoantigen prediction using `garnish_affinity`.
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
#' @seealso \code{\link{garnish_summary}}
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
#'\dontrun{
#'
#'  # load an example VCF
#'  dir <- system.file(package = "antigen.garnish") %>%
#'    file.path(., "extdata/testdata")
#'
#'    dt <- "antigen.garnish_example.vcf" %>%
#'    file.path(dir, .) %>%
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

garnish_variants <- function(vcfs, tumor_sample_name = "TUMOR"){

  # magrittr version check, this will not hide the error, only the NULL return on successful exit

  invisible(check_dep_versions())

  message("Loading VCFs")

sdt <- lapply(vcfs %>% seq_along, function(ivf){

  # load dt
      vcf <-  vcfR::read.vcfR(vcfs[ivf], verbose = TRUE)

        sample_id <- basename(vcfs[ivf])

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

  # rename generic columns to prevent downstream errors

  if (names(vdt) %like% "^V[0-9]+$" %>% any)
    vdt %>% data.table::setnames(names(vdt) %include% "^V[0-9]+$",
              paste(names(vdt) %include% "^V[0-9]+$", ".x", sep = ""))

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
    if (any(names(vdt) %like% "refseq")){

      vdt %<>% .[!cDNA_change %>% is.na &
                (!is.na(ensembl_transcript_id) | !is.na(refseq_id))]

      rs <- system.file(package = "antigen.garnish") %>% file.path(., "extdata", "Refseq_ids.txt") %>%
              data.table::fread

      vdt <- list(vdt[!is.na(ensembl_transcript_id)],
                  merge(vdt[is.na(ensembl_transcript_id)] %>% .[, ensembl_transcript_id := NULL],
                            rs, all.x = TRUE, by = "refseq_id")
              ) %>% data.table::rbindlist(use.names = TRUE)

      # drop redundancy from multiple NCBI tx ids matching to same ensembl tx id
      vdt %<>% .[!cDNA_change %>% is.na & !is.na(ensembl_transcript_id)] %>%
                  unique(by = c("sample_id", "cDNA_change", "ensembl_transcript_id"))

              }
      else{
        vdt %<>% .[!cDNA_change %>% is.na & !is.na(ensembl_transcript_id)]
      }

    if (vdt %>% nrow < 1){
      warning("No variants are present in the input file after filtering.")
      return(data.table::data.table(sample_id = sample_id))
      }

    if ("CF" %chin% (vdt %>% names))
      vdt %>% data.table::setnames("CF", "cellular_fraction")

    # check tumor_sample_name is in vcf
    # AF is a GT field not an INFO field, and account for multiple alt alleles in this scenario
    if (length(tumor_sample_name) != 1 | class(tumor_sample_name)[1] != "character")
      stop("A single tumor sample name must be provided.")

    afield <- paste(tumor_sample_name, "_AF", sep = "")

    if (!afield %chin% (vdt %>% names)){

      pns <- names(vdt) %include% "_AF" %>% stringr::str_replace("_AF",  "")


        # if sample names exist, report mismatches
      if (!(pns == "") %>% all){

        message("Unable to match tumor sample name, allelic_fraction will not be considered in downstream analysis.")
        message(paste("Possible sample names from the VCF are:", paste(pns, collapse = "\n"), sep = "\n"))
        message("To include allelic_fraction in downstream analysis, rerun garnish_variants with the correct tumor_sample_name argument.")
      }

    }

    if (afield %chin% (vdt %>% names)){

      vdt %>% data.table::setnames(afield, "allelic_fraction")

      vdt %<>% tidyr::separate_rows(c("ALT", "allelic_fraction"), sep = ",")

      # now keep only rows that match the previously split ANN field

      vdt <- vdt[ALT == stringr::str_extract(ANN, pattern = "^[AGCT]+(?=\\|)")]

    }

    return(vdt)

    })

  if(class(sdt)[1] == "list")
    sdt %<>% data.table::rbindlist(use.names = TRUE,
                                   fill = TRUE)

  if (nrow(sdt) == 0 || ncol(sdt) == 1){
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

  if ("cellular_fraction" %chin% names(sdt))
    sdt[, cellular_fraction := cellular_fraction %>% as.numeric]

  if ("allelic_fraction" %chin% names(sdt))
    sdt[, allelic_fraction := allelic_fraction %>% as.numeric]


  message("\nDone loading variants.")
  return(sdt)

}




#' Graph `garnish_affinity` results.
#'
#' Plot ADN, CDN, priority, frameshift, and fusion derived `nmers` for class I and class II MHC by `sample_id`.
#'
#' @param input Output from `garnish_affinity` to graph. `input` may be a data table object, list of data tables, or file path to a rio::import-compatible file type. If a list of data tables is provided, unique plots will be generated for each data table.
#' @param ext File extension passed to ggplot2::ggsave, default is "pdf".
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_summary}}
#'
#' @examples
#'\dontrun{
#'
#'# load example output
#'  dir <- system.file(package = "antigen.garnish") %>%
#'    file.path(., "extdata/testdata")
#'
#'    "antigen.garnish_PureCN_example_output.txt" %>%
#'    file.path(dir, .) %>%
#'
#'# create plot
#'    garnish_plot
#'
#'}
#'
#' @return NULL
#'
#' As a side effect: saves graph illustrating the number of neoantigens in each sample to the working directory. The threshold for inclusion of neoantigens is \eqn{IC_{50}} < 500nM.
#'
#' @export garnish_plot
#'
#' @md

garnish_plot <- function(input, ext = "pdf"){

  # magrittr version check, this will not hide the error, only the NULL return on successful exit
  invisible(check_dep_versions())

  ext <- paste0(".", ext)

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
              size = ggplot2::rel(2*1.425), face = "bold")) +
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

      class_col <- c("#ff6d00", "#2962ff", "#00c853", "#DE7AA7", "#BEBEBE")
      names(class_col) <- c("CDNs", "ADNs", "IEDB_high", "high_dissimilarity", "MHC_binders")

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

NA_to_0 <- function(v){

 if (!class(v)[1] %chin% c("numeric", "integer"))
    return(v)

    v[which(is.na(v))] <- 0

    return(v)

}

# summary plots

lapply(input %>% seq_along, function(i){

    dt <- input[[i]]

    # convert sample_id to character to prevent ggplot picking continuous X axis
    dt[, sample_id := as.character(sample_id)]

    s <- dt %>% garnish_summary

    n <- names(s)[which(names(s) %like% "CDNs|ADNs|MHC_binders|IEDB|dissimilarity")]

    s <- s[, .SD, .SDcols = c("sample_id", n)]

    s %<>% melt(id.vars = "sample_id")

    s[, class := variable %>% stringr::str_extract("class_(I$|II)")]

    s[class == "class_I", class := "MHC Class I"]
    s[class == "class_II", class := "MHC Class II"]

    s <- s[!is.na(class)] %>% .[!variable %like% "score"]

    s[, value := NA_to_0(value)]

    s[, variable := variable %>% stringr::str_replace("_class_.*$", "")]

    s %>% data.table::setnames(c("variable", "class", "value"),
                               c("type", "MHC", "N"))

  gdt <- s %>% gplot_missing_combn

  gg_dt <- data.table::rbindlist(list(s, gdt), use.names = TRUE)  %>%
        .[, N := max(N), by = c("sample_id", "type", "MHC")] %>% unique

  #check if only one class is present and if so drop that from faceting
  if (gg_dt[, all(N == 0), by  = "MHC"]$V1 %>% any){

    d <-gg_dt[, all(N == 0), by  = "MHC"]  %>%
              .[V1 == TRUE, MHC %>% unique]

    gg_dt <- gg_dt[MHC != d]

  }

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
            ggplot2::scale_fill_manual(values = class_col) +
            ggplot2::ggtitle(paste0("antigen.garnish summary")) +
            ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, keywidth = 1, keyheight = 1))

            ggplot2::ggsave(plot = g,
                  paste0("antigen.garnish_Neoepitopes_summary_",
                    gplot_fn,
                    ext)
                  , height = 6, width = 9)

    })

  return(NULL)

}



#' List top neoantigens for each sample and/or by clones within each sample.
#'
#' @param dt An output data table from `garnish_affinity`, either a data table object or path to a file.
#' @param nhits Integer. Maximum number of prioritized neoantigens per clone to return, default is 2.
#' @param binding_cutoff Numeric. Affinity threshold in nM to consider neoantigens, default is 500.
#'
#' @return A data table with the top neoantigen per sample, and if possible per clone, in rank order of clone frequency. Neoantigens are prioritized by permissive dissimilarity and iedb_score thresholds and then by differential agretopicity.
#'
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_summary}}
#'
#' @export garnish_antigens
#' @md

garnish_antigens <- function(dt, nhits = 2, binding_cutoff = 500){

  if (class(dt)[1] == "character") dt <- dt %>%
  rio::import %>%
    data.table::as.data.table

  if (class(dt)[1] == "data.frame") dt %<>%
    data.table::as.data.table

  dt %<>% data.table::copy

  if (!"Ensemble_score" %chin% names(dt))
    stop("Missing Ensemble_score column.  Input to garnish_antigens must be garnish_affinity output.")

  dt <- dt[Ensemble_score < binding_cutoff & pep_type != "wt"]

  c <- c("clone_id", "cl_proportion")

  if (!"clone_id" %chin% names(dt)){

    dt[, clone_id := as.numeric(NA)]

    c <- "clone_id"

  }

  if (!"min_DAI" %chin% names(dt) & "DAI" %chin% names(dt))
    dt[, min_DAI := DAI]

  ol <- c("dissimilarity", "iedb_score", "min_DAI", "Ensemble_score")

  ml <- ol[which(!ol %chin% names(dt))]

  ol <- ol[which(ol %chin% names(dt))]

  message(
    paste("Ranking peptides based on the follow metrics from available input: ",
    paste(ol, collapse = ", "), ".",
    sep = "")
  )

  if (length(ml) != 0) dt[, as.character(ml) := as.numeric(NA)]

  # split by  sample
  dtl <- dt %>% split(by = "sample_id")

  dtl <- lapply(dtl, function(n){

    dtll <- n %>% split(by = "clone_id")

    dtll <- lapply(dtll, function(nn){

      if (nrow(nn[dissimilarity > 0 & iedb_score > 0]) > 0)
        return(nn[order(min_DAI, decreasing = TRUE)][1:nhits])

      if (nrow(nn[dissimilarity > 0]) > 0)
        return(nn[order(min_DAI, decreasing = TRUE)][1:nhits])

      if (nrow(nn[iedb_score > 0]) > 0)
        return(nn[order(min_DAI, decreasing = TRUE)][1:nhits])

      if (nrow(nn[min_DAI > 1]) > 0)
          return(nn[order(min_DAI, decreasing = TRUE)][1:nhits])

      return(nn[order(Ensemble_score, decreasing = FALSE)][1:nhits])

    }) %>% data.table::rbindlist(use.names = TRUE)

  }) %>% data.table::rbindlist(use.names = TRUE)

  # now drop extra rows added by nhits argument which will all have NA values
  dtl <- dtl[!is.na(Ensemble_score)]

  n <- names(dtl)[which(names(dtl) %chin% c("cDNA_change", "protein_change", "external_gene_name"))]

  if (length(n) < 1) n <- NULL

  dt <- dtl[, .SD %>% unique, .SDcols = c("sample_id", "nmer", "MHC", n,
                                        "Ensemble_score", "dissimilarity", "iedb_score", "min_DAI", c)]

  if (!"clone_id" %chin% names(dt)) dt <- dt %>% .[order(sample_id)]

  if ("clone_id" %chin% names(dt)) dt <- dt %>% .[order(sample_id, clone_id)]

  return(dt)

}
