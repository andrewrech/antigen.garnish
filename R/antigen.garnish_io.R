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
#' @return A summary data table of \code{dt} by \code{sample_id} with the following columns:
#' * **priority_neos**: neoepitopes that meet both classic and alternative criteria, or that meet classic criteria and are derived from frameshift mutations
#' * **classic_neos**: high affinity neoepitopes
#' * **classic_top_score**: sum of the top three affinity scores (1 / MHC affinity (nM))
#' * **alt_neos**: mutant nmers predicted to bind MHC with greatly improved affinity relative to their non-mutated counterpart
#' * **alt_top_score**: sum of the top three `nmer` differential agretopicity indices; differential agretopicity index (DAI) is the ratio of MHC binding affinity between mutant and wt peptide.
#' * **nmers**: mutant nmers
#' * **variants**: genomic variants evaluated
#' * **transcripts**: transcripts evaluated
#' * **predictions**: wt and mutant predictions performed
#' * **mhc_binders**: nmers predicted to at least minimally bind MHC
#' @export garnish_summary
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
                          mutant_loc > (mismatch_s + 2)
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
garnish_variants <- function(vcfs) {

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

    if(vdt %>% nrow < 1) return(data.table::data.table(sample_id = sample_id))

    # filter passing Strelka2 and muTect variants
    if(vcf_type == "Strelka") vdt <- vdt[FILTER == "PASS"]
    if(vcf_type == "Mutect") vdt <- vdt[INFO %>%
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
