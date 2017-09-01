## ---- garnish_summary
#' Summarize epitope prediction.
#'
#' Calculate neoepitope summary statistics over samples.
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
#'    garnish_variants %>%
#'
#'  # predict neoepitopes
#'    garnish_predictions %>%
#'
#'  # summarize predictions
#'    garnish_summary %T>%
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
#' * **nmers**: total mutant `nmers` created
#' * **variants**: total genomic variants evaluated
#' * **transcripts**: total transcripts evaluated
#' * **predictions**: wt and mutant predictions performed
#' * **mhc_binders**: `nmers` predicted to at least minimally bind MHC (< 5000nM \eqn{IC_{50}})
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

  dt <- dt[DAI != Inf & DAI != -Inf & Consensus_scores != Inf & Consensus_scores != -Inf]

    # function to sum top values of a numeric vector

  sum_top_v <- function(x, value = 3){

        x %<>% sort %>% rev
        return(sum(x[1:value]))
      }

  dt %>% data.table::setkey(sample_id)

  dtn <- parallel::mclapply(dt[, sample_id %>% unique], function(id){

    dt <- dt[sample_id == id]

      return(
          data.table::data.table(
          sample_id = id,
          priority_neos = dt[DAI > 10 & (
                       Consensus_scores < 50 |
                        (frameshift == TRUE &
                         (
                          mutant_index > (mismatch_s + 2)
                          )
                          )
                       )] %>% nrow,
          classic_neos = dt[Consensus_scores < 50] %>% nrow,
          classic_top_score = dt[Consensus_scores < 5000, (1/Consensus_scores) %>% sum_top_v],
          alt_neos = dt[Consensus_scores < 5000 & DAI > 10] %>% nrow,
          alt_top_score = dt[Consensus_scores < 5000, DAI %>% sum_top_v],
          mhc_binders = dt[Consensus_scores < 5000] %>% nrow,
          variants = dt[, var_uuid %>% unique] %>% length,
          transcripts = dt[ensembl_transcript_id %>% unique] %>% nrow,
          predictions = dt[pep_type %like% "mut"] %>% nrow,
          nmers = dt[pep_type %like% "mut", nmer %>% unique] %>% length
          ))

    }) %>% data.table::rbindlist

  return(dtn)
}



## ---- garnish_variants
#' Intakes variants and returns an intersected data table for epitope prediction.
#'
#' Process variants annotated with [SnpEff](http://snpeff.sourceforge.net/). VCFs from matched samples are optionally intersected. [MuTect2](https://github.com/broadinstitute/gatk)/[Strelka2](https://github.com/Illumina/strelka)-derived VCFs are filtered for high confidence variants prior to intersection.
#'
#' @param vcfs Character vector. VFC files to import.
#'
#' @return A VCF as a data table with one unique SnpEff variant annotation per row, including:
#' * **sample_id**: sample identifier constructed from input \code{.bam} file names
#' * **se**: SnpEff annotation
#' * **effect_type**: SnpEff effect type
#' * **ensembl_transcript_id**: transcript effected
#' * **ensembl_gene_id**: gene effected
#' * **protein_change**: protein change in [HGVS](http://varnomen.hgvs.org/recommendations/DNA/) format
#' * **cDNA_change**: cDNA change in [HGVS](http://varnomen.hgvs.org/recommendations/protein/) format
#' * **protein_coding**: is the variant protein coding?
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

    vdt %<>% get_snpeff_annot

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
