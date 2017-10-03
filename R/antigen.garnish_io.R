## ---- garnish_summary
#' Summarize neoepitope prediction.
#'
#' Calculate neoepitope summary statistics for priority, classic, and alternative neoepitopes for each sample.
#'
#' @param dt Data table. Prediction output from `garnish_predictions`.
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
#'    antigen.garnish::garnish_predictions %>%
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
#' * **alt_top_score**: sum of the top three mutant `nmer` differential agretopicity indices; differential agretopicity index (DAI) is the ratio of MHC binding affinity between mutant and wt peptide (see below)
#' * **priority_neos**: mutant peptides that meet both ADN and CDN criteria, or that meet CDN criteria and are derived from frameshift mutations
#' * **fs_neos**: mutant `nmers` derived from frameshift variants predicted to bind MHC with < 500nM \eqn{IC_{50}})
#' * **fusion_neos**: mutant `nmers` derived from fusion variants predicted to bind MHC with < 500nM \eqn{IC_{50}})
#' * **nmers**: total mutant `nmers` created
#' * **predictions**: wt and mutant predictions performed
#' * **mhc_binders**: `nmers` predicted to at least minimally bind MHC (< 5000nM \eqn{IC_{50}})
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
#'To better model potential for oligoclonal antitumor responses directed against neoepitopes, we additionally report a **top three neoepitope score**, which is defined as the sum of the top three affinity scores \eqn{\left(\frac{1}{IC_{50}}\right)} for CDNs or sum of top three DAI for ADNs. The top three was chosen in each case because this is the minimum number that captures the potential for an oligoclonal T cell response and mirrors experimentally confirmed oligoclonality of T cell responses against human tumors. Moreover, the top three score was the least correlated to total neoepitope load (vs. top 4 through top 15) is a large scale human analysis of neoepitope across 27 disease types (R-squared = 0.0495), and therefore not purely a derivative of total neoepitope load.
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

# summarize over unique nmers

  dt %<>% data.table::as.data.table %>%
    data.table::copy %>%
    unique(by = c("pep_type",
                  "MHC",
                  "nmer"))

  # set NA DAI to 0 to filter Inf and -Inf
    dt[, DAI := DAI %>% as.numeric]
    dt[is.na(DAI), DAI  := 0]

  dt <- data.table(DAI = NA)
  dt <- dt[DAI != Inf]

  dt <- dt[DAI != Inf & DAI != -Inf & Consensus_scores != Inf & Consensus_scores != -Inf]

# function to sum top values of a numeric vector

  sum_top_v <- function(x, value = 3){

        x %<>% sort %>% rev
        return(sum(x[1:value]))
      }

  dt %>% data.table::setkey(sample_id)

  dtn <- parallel::mclapply(dt[, sample_id %>% unique], function(id){

    dt <- dt[sample_id == id]

    if (!("fusion_uuid" %chin% names(dt)) %>% any) dt[, fusion_uuid := NA]
    if (!("effect_type" %chin% names(dt)) %>% any) dt[, effect_type := NA]

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
                                      Consensus_scores < 5000, (1/Consensus_scores) %>% sum_top_v],
          classic_top_score_class_II = dt[class == "II" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 5000, (1/Consensus_scores) %>% sum_top_v],
          alt_neos_class_I = dt[class == "I" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 5000 &
                                      DAI > 10] %>%
                                      data.table::as.data.table %>% nrow,
          alt_neos_class_II = dt[class == "II" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 5000 &
                                      DAI > 10] %>%
                                      data.table::as.data.table %>%  nrow,
          alt_top_score_class_I = dt[class == "I" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 5000, DAI %>% sum_top_v],
          alt_top_score_class_II = dt[class == "II" &
                                      pep_type == "mutnfs" &
                                      Consensus_scores < 5000, DAI %>% sum_top_v],
          fs_neos_class_I = dt[class == "I" &
                                      effect_type %like% "frameshift" &
                                      Consensus_scores < 500] %>%
                                      data.table::as.data.table %>% nrow,
          fs_neos_class_II = dt[class == "II" &
                                      effect_type %like% "frameshift" &
                                      Consensus_scores < 500]  %>%
                                      data.table::as.data.table %>% nrow,
          fusion_neos_class_I = dt[class == "I" &
                                      !is.na(fusion_uuid) &
                                      Consensus_scores < 500] %>%
                                      data.table::as.data.table %>% nrow,
          fusion_neos_class_II = dt[class == "II" &
                                      !is.na(fusion_uuid) &
                                      Consensus_scores < 500]  %>%
                                      data.table::as.data.table %>% nrow,
          mhc_binders_class_I = dt[class == "I" &
                                      Consensus_scores < 5000] %>% nrow,
          mhc_binders_class_II = dt[class == "II" &
                                      Consensus_scores < 5000] %>% nrow,
          predictions = dt[pep_type %like% "mut"] %>% nrow,
          nmers = dt[pep_type %like% "mut", nmer %>% unique] %>% length
          ))

    }) %>% data.table::rbindlist

# get variant level statistics if available
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
#' Process VCF variants annotated with [SnpEff](http://snpeff.sourceforge.net/) for neoepitope prediction using `garnish_predictions`. VCFs from matched samples are optionally intersected. [MuTect2](https://github.com/broadinstitute/gatk)/[Strelka2](https://github.com/Illumina/strelka)-derived VCFs are filtered for high confidence variants prior to intersection.
#'
#' @param vcfs Paths to VFC files to import.
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
#' @seealso \code{\link{garnish_predictions}}
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
#' @export garnish_variants
#' @md

garnish_variants <- function(vcfs){

  message("Loading VCFs")

  ivfdtl <- parallel::mclapply(vcfs %>% seq_along, function(ivf){

  # load dt
      vcf <-  vcfR::read.vcfR(vcfs[ivf], verbose = TRUE)

  # extract sample names from Mutect2 and Strelka command line for intersection

      sample_id <- vcf@meta %>%
                      stringr::str_extract_all("[^ ]+\\.bam") %>%
                      unlist %>%
                      unique %>%
                      basename %>%
                      sort %>%
                      paste(collapse = "_") %>%
                      stringr::str_replace("\\.bam", "")
      # extract vcf type
      vcf_type <- vcf@meta %>% unlist %>% stringr::str_extract(stringr::regex("Strelka|Mutect", ignore_case = TRUE)) %>% stats::na.omit %>%
                  unlist %>% data.table::first
      if (vcf_type %>% length == 0) vcf_type <- "other"

  # return a data table of variants

    vdt <- vcf@fix %>% data.table::as.data.table

    if (vcf@gt %>% length > 0) vdt <- cbind(vdt, vcf@gt %>% data.table::as.data.table)

    if (vdt %>% nrow < 1) return(data.table::data.table(sample_id = sample_id))

    # filter passing Strelka2 and muTect variants
    if (vcf_type == "Strelka") vdt <- vdt[FILTER == "PASS"]
    if (vcf_type == "Mutect") vdt <- vdt[INFO %>%
                                        stringr::str_extract("(?<=TLOD=)[0-9\\.]") %>%
                                        as.numeric > 6.0]
    vdt[, sample_id := sample_id]
    vdt[, vcf_type := vcf_type]

    # filter SnpEff warnings

    vdt %<>% get_snpeff

    vdt %<>% .[!se %like% "ERROR_.*CHROMOSOME"]
    vdt %<>% .[!se %likef% "WARNING_SEQUENCE_NOT_AVAILABLE"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_INCOMPLETE"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_NO_START_CODON"]

    return(vdt)
    })

    ivfdt <- ivfdtl %>% data.table::rbindlist

  merge_vcf <- function(dt, dt2){

    # a function to intersect annotated variants across VCFs using SnpEff

      sdt <- merge(dt, dt2[, .SD,
             .SDcols = c("CHROM",
              "POS", "REF", "cDNA_change")],
               by = c("CHROM", "POS", "REF",
              "cDNA_change")) %>% .[, vcf_type := "intersect"]
      return(sdt)

    }
  # return an intersected data table of variants

  sdt <- parallel::mclapply(ivfdt[, sample_id %>% unique], function(sn){

    # find data tables with matching sample names

  sdt <- lapply(ivfdtl, function(dt){

     dt[, sample_id %>% .[1]] == sn

    }) %>% unlist

  # merge all data tables with matching sample names
    if (ivfdtl[sdt] %>% length == 1) return(ivfdtl[[sdt %>% which]])
    if (ivfdtl[sdt] %>% length > 1) return(ivfdtl[sdt] %>% Reduce(merge_vcf, .))


  }) %>% data.table::rbindlist

  # select protein coding variants without NA
  sdt %<>% .[protein_coding == TRUE &
            !protein_change %>% is.na &
            !effect_type %>% is.na &
            effect_type %like% "insertion|deletion|missense|frameshift"]

  return(sdt)

}



## ---- garnish_plot
#' Graph `garnish_summary` results.
#'
#' Plot ADN, CDN, priority, frameshift, and fusion derived `nmers` for class I and class II MHC by `sample_id`.
#'
#' @param dt Data table. Non gene-fusion prediction output from `garnish_predictions`. Either dt or jdt must be provided as input.
#' @param jdt Data table. Prediction output from `garnish_jaffa` -> `garnish_predictions`. Either dt or jdt must be provided as input.
#' @param save Optional: Should the plot be saved to the working directory as antigen.garnish_summary.pdf? Default is TRUE.
#'
#' @seealso \code{\link{garnish_predictions}}
#' @seealso \code{\link{garnish_summary}}
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download an example VCF
#'    g <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'  # extract variants
#'    antigen.garnish::garnish_variants %>%
#'
#'  # add test MHC types
#'      .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
#'                   "H-2-Kb H-2-IAd",
#'                  "HLA-A*01:47 HLA-DRB1*03:07")] %>%
#'
#'  # predict neoepitopes
#'    antigen.garnish::garnish_predictions %>%
#'
#'  # plot it
#'    antigen.garnish::garnish_plot
#'
#'  # here is our ggplot object that can be further customized
#'   class(g)
#'
#'  # save the plot with a custom ggplot theme
#'  g + my_ggplot_custom_theme_object %>% cowplot::ggsave("antigen.garnish_custom_themed_plot.pdf")
#'
#'  # show location of where default plot was saved
#'   list.files(pattern = "summary\\.pdf", full.names = TRUE)
#'
#'}
#'
#' @return
#'
#' A list of ggplot2 objects and, if save = TRUE, saves each of these objects pdfs in the working directory.
#' The first plot is a column graph indicating the number of peptides in each sample classified as ADN, CDN, priority, frameshift-derived, and fusion-derived neoepitopes
#' (if any) for each sample, faceted by Class I vs Class II MHC. The threshold for fusion and frameshift derived neoepitopes is Consensus_scores < 1000nM.
#' See `garnish_summary`` documentation for further explanation of neoepitope classification of ADN, CDN, and priority.
#'
#' The second plot and third plots are returned only if frameshift and/or fusion mutants (respectively) are present in the input table. These plot
#' the number of peptides per sample by MHC class I or II and binned by binding affinity (1000 - 500nM, 500 - 50nM, and < 50 nM).
#'
#' @export garnish_plot
#'
#' @md

garnish_plot <- function(dt = NULL, jdt = NULL, save = TRUE){

  ag_gg_theme <-
    ggplot2::theme(line = ggplot2::element_line(colour = "#000000")) +
    ggplot2::theme(axis.line = ggplot2::element_line(color="black")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(2*1.2),
          colour = "#000000", lineheight = 0.9, face = "bold", vjust = 0)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "#000000",
          size = ggplot2::rel(2*1.425), face = "bold")) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(colour = "#000000",
          size = ggplot2::rel(2*1.425), face = "bold", angle = 90, vjust = 0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(colour = "#000000",
          size = ggplot2::rel(2*1.425), face = "bold", angle = 30, hjust = 0.9, vjust = 0.92)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(colour = "#000000",
          size = ggplot2::rel(2*1.425), face = "bold", hjust = 0.9, vjust = 0.92, angle = 30)) +
    ggplot2::theme(legend.text = ggplot2::element_text(colour = "#000000",
          size = ggplot2::rel(2*1))) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "#eeeeee")) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold")) +
    ggplot2::theme(legend.key = ggplot2::element_blank()) +
    ggplot2::theme(legend.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
    ggplot2::theme(plot.margin=grid::unit(c(0, 0, 0, 0 ), "cm"))

  ag_colors <-   c("#ff80ab",
                   "#b388ff",
                   "#82b1ff",
                   "#a7ffeb",
                   "#b9f6ca",
                   "#f4ff81",
                   "#ffe57f",
                   "#ff9e80",
                   "#ff8a80",
                   "#ea80fc",
                   "#8c9eff",
                   "#80d8ff",
                   "#84ffff",
                   "#ccff90",
                   "#ffff8d",
                   "#ffd180")

  # check function input
    if (missing(dt) & missing(jdt)) stop("At least one dt from garnish_predictions must be provided.")

    if (!missing(dt) & !missing(jdt)) dt <- cat_tables(dt, jdt)

    if(!missing(dt)){

      if (!(c("nmer",
              "MHC",
              "sample_id",
              "DAI",
              "Consensus_scores") %chin% names(dt)) %>% any)
        stop("'sample_id', 'nmer', 'MHC', 'frameshift', 'DAI', and 'Consensus_scores' columns are required in dt")
  # create and filter data table
    dt <- dt[pep_type != "wt"] %>% unique(by = c("nmer",
                                                 "MHC",
                                                 "sample_id",
                                                 "DAI",
                                                 "Consensus_scores"))

    dt <- dt %>% .[Consensus_scores < 5000] %>%
      .[pep_type == "mutnfs" & Consensus_scores < 50, type := "CDN"] %>%
      .[pep_type == "mutnfs" & DAI > 10, type := "ADN"] %>%
      .[pep_type == "mutnfs" & Consensus_scores < 50 & DAI > 10, type := "priority"] %>%
      .[pep_type == "fus" & Consensus_scores < 1000, type := "fusion"] %>%
      .[pep_type == "mut_other" & Consensus_scores < 1000, type := "frameshift"]

    if (nrow(dt) < 1)
      stop("No neoeptiopes with Consensus_scores < 5000nM ")
    if (nrow(dt_pl) < 1)
      stop("No neoepitopes meet classification criteria")

    gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "type")]

  # function to graph dt

    gdt <- dt_pl %>% (function(dt){

      ns <- dt[, sample_id %>% unique %>% length]

      gdt <- data.table(sample_id = dt[, sample_id %>% unique],
                        MHC = c(replicate(ns, "Class I"), replicate(ns, "Class II")),
                        type = c(replicate(ns * 2, "ADN"), replicate(ns * 2, "CDN"), replicate(ns * 2, "priority"),
                                 replicate(ns * 2, "frameshift"), replicate(ns * 2, "fusion")),
                        N = 0) %>% unique
      return(gdt)
    })

    gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>%
      .[, N := max(N), by = c("sample_id", "type", "MHC")] %>% unique

  # shorten names for display

    if (any(gg_dt[, sample_id %>% unique %>% nchar] > 7)){

      for (i in gg_dt[nchar(sample_id) > 7, sample_id %>% unique] %>% seq_along){

        gg_dt[sample_id == gg_dt[nchar(sample_id) > 7, sample_id %>% unique][i],
              sample_id := sample_id %>%
                substr(1, 7) %>% paste0(., "_", i)]
      }
    }

    g1 <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
            ggplot2::geom_col(ggplot2::aes(fill = type), col = "black", position = "dodge") +
            ggplot2::facet_wrap(~MHC) +
            ag_gg_theme +
            ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
            ggplot2::scale_fill_manual(values = ag_colors) +
            ggplot2::theme(strip.text.x = ggplot2::element_text(size = ggplot2::rel(2))) +
            ggplot2::ylab("peptides") +
            ggplot2::xlab("") +
            ggplot2::ggtitle(paste0("antigen.garnish summary"))

    cowplot::ggsave(plot = g1, "antigen.garnish_summary.pdf", height = 6, width = 9)

    if (!(dt[, effect_type %>% unique] %like% "frameshift" %>% any)) g2 <- NA

  # frameshift plot

    if (dt[, effect_type %>% unique] %like% "frameshift" %>% any){

      dt_pl <- dt[MHC %like% "(HLA-[ABC]\\*)|(H-2-[A-Z][a-z])", MHC := "Class I"] %>%
        .[MHC %like% "(HLA-D[A-Z0-9]+\\*)|(H-2-[A-Z]{2}[a-z])", MHC := "Class II"] %>%
        .[effect_type %like% "frameshift" & Consensus_scores < 1000]

      dt_pl[, binding := "<1000nM"] %>%
        .[Consensus_scores < 500, binding := "<500nM"] %>%
        .[Consensus_scores < 50, binding := "<50nM"]

      gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "binding")]

      gdt <- dt_pl %>% (function(dt){

        ns <- dt[, sample_id %>% unique %>% length]

        gdt <- data.table(sample_id = dt[, sample_id %>% unique],
                          MHC = c(replicate(ns, "Class I"), replicate(ns, "Class II")),
                          binding = c(replicate(ns * 2, "<50nM"), replicate(ns * 2, "<500nM"), replicate(ns * 2, "<1000nM")),
                          N = 0) %>% unique
      })

      gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>%
        .[, N := max(N), by = c("sample_id", "binding", "MHC")] %>% unique

      if (any(gg_dt[, sample_id %>% unique %>% nchar] > 7)){

  # shorten names for display

        for (i in gg_dt[nchar(sample_id) > 7, sample_id %>% unique] %>% seq_along){

          gg_dt[sample_id == gg_dt[nchar(sample_id) > 7, sample_id %>% unique][i],
                sample_id := sample_id %>%
                  substr(1, 7) %>% paste0(., "_", i)]
        }
      }

      g2 <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
              ggplot2::geom_col(ggplot2::aes(fill = binding), col = "black", position = "dodge") +
              ggplot2::facet_wrap(~MHC) +
              ag_gg_theme +
              ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
              ggplot2::scale_fill_manual(values = ag_colors[1:3]) +
              ggplot2::theme(strip.text.x = ggplot2::element_text(size = ggplot2::rel(2))) +
              ggplot2::ylab("peptides") +
              ggplot2::xlab("") +
              ggplot2::ggtitle(paste0("Frameshift neoepitopes"))

    cowplot::ggsave(plot = g2, "antigen.garnish_Frameshifts_summary.pdf", height = 6, width = 9)

    }
}

  # fusion plot

    if(!missing(jdt) & missing(dt)) jdt <- jdt

    if(!missing(jdt) & !missing(dt)) jdt <- dt

    if (!"fusion genes" %chin% (jdt %>% names)) g3 <- NA

    if ("fusion genes" %chin% (jdt %>% names)){


    dt_pl <- jdt[MHC %like% "(HLA-[ABC]\\*)|(H-2-[A-Z][a-z])", MHC := "Class I"] %>%
      .[MHC %like% "(HLA-D[A-Z0-9]+\\*)|(H-2-[A-Z]{2}[a-z])", MHC := "Class II"] %>%
      .[!is.na(fusion_uuid) & Consensus_scores < 1000]

    dt_pl[, binding := "<1000nM"] %>%
      .[Consensus_scores < 500, binding := "<500nM"] %>%
      .[Consensus_scores < 50, binding := "<50nM"]

    gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "binding")]

    gdt <- dt_pl %>% (function(dt){

      ns <- dt[, sample_id %>% unique %>% length]

      gdt <- data.table(sample_id = dt[, sample_id %>% unique],
                        MHC = c(replicate(ns, "Class I"), replicate(ns, "Class II")),
                        binding = c(replicate(ns * 2, "<50nM"), replicate(ns * 2, "<500nM"), replicate(ns * 2, "<1000nM")),
                        N = 0) %>% unique
    })

    gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>%
      .[, N := max(N), by = c("sample_id", "binding", "MHC")] %>% unique

    if (any(gg_dt[, sample_id %>% unique %>% nchar] > 7)){

      message("Sample_id has >7 characters, shortening names for aesthetics.
              To circumvent this, change sample_ids to less than 7 characters in the input data.table.")

      for (i in gg_dt[nchar(sample_id) > 7, sample_id %>% unique] %>% seq_along){

        gg_dt[sample_id == gg_dt[nchar(sample_id) > 7, sample_id %>% unique][i],
              sample_id := sample_id %>%
                substr(1, 7) %>% paste0(., "_", i)]
      }
    }

    g3 <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
      ggplot2::geom_col(ggplot2::aes(fill = binding), col = "black", position = "dodge") +
      ggplot2::facet_wrap(~MHC) +
      ag_gg_theme +
      ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
      ggplot2::scale_fill_manual(values = ag_colors[1:3]) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = ggplot2::rel(2))) +
      ggplot2::ylab("peptides") +
      ggplot2::xlab("") +
      ggplot2::ggtitle(paste0("Fusion neoepitopes"))

   cowplot::ggsave(plot = g3, "antigen.garnish_Fusions_summary.pdf", height = 6, width = 9)

  }

  if (missing(dt) & !missing(jdt)) return(g3)

  if (missing(jdt) & !missing(dt)) return(list(g1,g2))

  glist <- list(g1, g2, g3)

  return(glist)

  }



## ---- cat_tables
#' Internal function to combine `garnish_predictions` outputs from garnish_jaffa
#' and garnish_variants or direct input to garnish_predictions.
#'
#' @param vdt WRITE DESC
#' @param jdt WRITE DESC
#' @export cat_tables
#' @md
#'
cat_tables <- function(vdt, jdt){

  # get the jaffa_dt

    ifelse ("fus_tx" %chin% names(dt1), {
      jdt <- dt1

    })
    if ("fus_tx" %chin% names(fus_tx) %>% any) jdt <- dt1
    if ("vdt" %chin% names(dt2) %>% any) jdt <- dt2

  # get the vcf predictions output dt

    if (grepl("mutnfs", (dt1[, .SD %>% unique, .SDcols = "pep_type"] %>%
                         unlist)) %>% any) vdt <- dt1
    if (grepl("mutnfs", (dt2[, .SD %>% unique, .SDcols = "pep_type"] %>%
                         unlist)) %>% any) vdt <- dt2

  # subset tables for merge

    vdt <- vdt[, .SD, .SDcols = c("nmer", "sample_id", "MHC", "Consensus_scores", "DAI",
                                  "Lower.CI", "Upper.CI", "ensembl_transcript_id",
                                  "effect_type", "external_gene_name", "frameshift", "protein_change",
                                  "cDNA_change", "mutant_index", "pep_mut", "pep_wt", "pep_type", "var_uuid", "dai_uuid", "nmer_uuid")]

    jdt <- jdt[, .SD, .SDcols = c("nmer", "sample_id", "MHC", "Consensus_scores", "Lower.CI", "Upper.CI",
                                  "fusion genes", "frameshift", "mutant_index", "pep_mut", "pep_gene_1", "var_uuid",
                                  "fusion_uuid", "nmer_uuid")]

  # add fus as pep_type for later analysis, add mutfs for frameshifts
    jdt[, pep_type := "fus"]

    vdt[is.na(pep_type) & effect_type %like% "frameshift",
        pep_type := "mutfs"]

  # change names to overlap in the vcf dt
    strings <- jdt[, sample_id %>% unique]

  # contruct our regex
   strings <- paste("(", strings, ")", sep = "", collapse = "|")

  # rename stings in vcf dt from full BAM name
    vdt[, sample_id := sample_id %>% stringr::str_extract(pattern = strings)]

  # final merge
    dto <- merge(vdt, jdt, all = TRUE, by = intersect(names(vdt), names(jdt)))
    return(dto)
}

