## ---- garnish_summary
#' Summarize neoepitope prediction.
#'
#' Calculate neoepitope summary statistics for priority, classic, alternative, frameshift-derived, and fusion-derived neoepitopes for each sample.
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
#' * **alt_top_score**: sum of the top three mutant `nmer` differential agretopicity indices; differential agretopicity index (DAI) is the ratio of MHC binding affinity between mutant and wt peptide (see below). If the peptide is fusion- or frameshift-derived, the binding affinity of the closest wt peptide determined by BLAST is used to calculate DAI.
#' * **priority_neos**: mutant peptides that meet both ADN and CDN criteria, and if applicable, have a fitness_score of >= 1.
#' * **fs_neos**: mutant `nmers` derived from frameshift variants predicted to bind MHC with < 500nM \eqn{IC_{50}}
#' * **fusion_neos**: mutant `nmers` derived from fusion variants predicted to bind MHC with < 500nM \eqn{IC_{50}}
#' * **nmers**: total mutant `nmers` created
#' * **predictions**: wt and mutant predictions performed
#' * **mhc_binders**: `nmers` predicted to at least minimally bind MHC (< 5000nM \eqn{IC_{50}})
#' * **fitness_scores**: Sum of the top 3 fitness_score values per sample. See `?garnish_fitness`.
#' * **garnish_score**: Sample level immune fitness parameter modeled by summing the exponential of fitness_scores for each dominant neoepitope across all clones in the tumor sample. See ?garnish_predicitons.
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
#' @seealso \code{\link{garnish_plot}} for a sample level summary plotting function.
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


  dt <- dt[DAI != Inf & DAI != -Inf & Consensus_scores != Inf & Consensus_scores != -Inf]

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
    dtn <- merge(dtn, dt[, .SD %>% unique, .SDcols = c("sample_id", "garnish_score")],
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
#' Process paired tumor-normal VCF variants annotated with [SnpEff](http://snpeff.sourceforge.net/) for neoepitope prediction using `garnish_predictions`. VCFs from matched samples can be optionally intersected to select only variants present across input files.
#'
#' Recommended somatic variant callers: [MuTect2](https://github.com/broadinstitute/gatk), [Strelka2](https://github.com/Illumina/strelka)
#'
#' **Multi-sample vcfs are not supported**.  Single sample and paired tumor-normal vcfs are required.  Multi-sample vcfs will be considered a single sample in the output.
#' The prop_tab argument may be used to provide clonality/cell fraction information with variants, or allelic fractions.  If the vcf is annotated in the info field and this is indicated in the header, the value "CF" will be used to obtain the cell fraction for variants.
#'
#' @param vcfs Paths to one or more VFC files to import.
#' @param intersect Logical. Return only the intersection of variants in multiple `vcfs` with identical sample names? Intersection performed on `SnpEff` annotations. One `vcf` file per somatic variant caller-input samples pair is required.
#' @param prop_tab. Optional.  Character vector of file paths to a table with clonality or allelic frequencies in the same order as corresponding VCFs. If allelic frequencies are given clonality will be inferred with the assumptions of even ploidy, monophylogeny, no loss of truncal mutations, and all allele frequencies lower than the maximum per sample represent subclonal mutations, this method has not been validated. Acceptable formats for input include:
#'
#'dt with clonality:
#'
#'
#'     Column name                 Example input
#'
#'     chr                         X
#'     start                       4550159
#'     end                         4550159
#'     CELLFRACTION                0.15
#'
#'dt with allele fractions:
#'
#'
#'     Column name                 Example input
#'
#'     #CHROM                      X
#'     POS                         4550159
#'     REF                         A
#'     ALT                         GC
#'     AF                          0.85
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
#' if prop_tab is provided either:
#' * **CELLFRACTION**: cell fraction taken from input, such as from clonality estimates from [PureCN](http://www.github.com/lima1/PureCN)
#' * **AF**: allelic fraction taken from input
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
#'
#' @references
#' Callari M, Sammut SJ, De Mattos-Arruda L, Bruna A, Rueda OM, Chin SF, and Caldas C. 2017. Intersect-then-combine approach: improving the performance of somatic variant calling in whole exome sequencing data using multiple aligners and callers. Genome Medicine. 9:35.
#'
#' @export garnish_variants
#' @md

garnish_variants <- function(vcfs, intersect = TRUE, prop_tab = NULL){

  # magrittr version check, this will not hide the error, only the NULL return on successful exit
  invisible(check_dep_versions())

  message("Loading VCFs")

  if (!missing(prop_tab) & intersect == TRUE){

    message("prop_tab argument is present, setting `intersect` to FALSE.
    To intersect and then use clonality, add an \"AF\" or \"CELLFRACTION\" to the table before passing to `garnish_predictions`.")
    intersect <- FALSE

  }

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

  # account for no bam names in header, as with VarScan
      if (sample_id == "") sample_id <- basename(vcfs[ivf])

  # extract vcf type
      vcf_type <- vcf@meta %>%
        unlist %>%
        stringr::str_extract(stringr::regex("(Strelka)|(Mutect)|(VarScan)|(samtools mpileup)|(somaticsniper)|(freebayes)|(virmid)",
          ignore_case = TRUE)) %>%
        stats::na.omit %>%
                 unlist %>%
        data.table::first
      if (vcf_type %>% length == 0) vcf_type <- "other"

  # return a data table of variants

    vdt <- vcf@fix %>% data.table::as.data.table

  # check that VCF is SnpEff-annotated

    if (
        vdt[, INFO %>% unique] %>% is.na ||
        !vdt[, INFO %>%
      stringr::str_detect(stringr::fixed("ANN=")) %>% all]
      )
      stop(paste0("\nInput file \n", vcfs[ivf], "\nis missing INFO SnpEff annotations"))

    if (vcf@gt %>% length > 0) vdt <- cbind(vdt, vcf@gt %>% data.table::as.data.table)

    if (vdt %>% nrow < 1) return(data.table::data.table(sample_id = sample_id))

    vdt[, sample_id := sample_id]
    vdt[, vcf_type := vcf_type]

    # filter SnpEff warnings

    vdt %<>% get_snpeff

    vdt %<>% .[!se %like% "ERROR_.*CHROMOSOME"]
    vdt %<>% .[!se %likef% "WARNING_SEQUENCE_NOT_AVAILABLE"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_INCOMPLETE"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_NO_START_CODON"]

    if (vdt %>% nrow < 1) return(data.table::data.table(sample_id = sample_id))

    # filter out NA
    vdt %<>% .[!ensembl_transcript_id %>% is.na &
               !cDNA_change %>% is.na]

    # this bugs downstream if nrow = 0 at this point, ie vcf of all intergenic
    if (vdt %>% nrow < 1) return(data.table::data.table(sample_id = sample_id))

    if (any(stringr::str_detect(vcf@meta, stringr::fixed("ID=CF"))))
      vdt[, CELLFRACTION := INFO %>% stringr::str_extract("(?<=(CF\\=))[01]\\.[0-9]+(?=;)") %>% as.numeric]

    return(vdt)
    })

rename_vcf_fields <- function(dt1, dt2) {

       # a function to rename VCF INFO and
        # FORMAT fields for merging

      for (common_col in c("INFO", "FORMAT")){

          if (common_col %chin% (dt %>% names))
            dt %>%
              data.table::setnames(common_col,
                paste0(common_col, "_", dt[, vcf_type[1]]))
          if (common_col %chin% (dt2 %>% names))
            dt2 %>%
              data.table::setnames(common_col,
                paste0(common_col, "_", dt2[, vcf_type[1]]))
        }
    }

merge_vcf <- function(dt, dt2){

    # a function to intersect annotated variants
      # across VCFs using SnpEff

      sdt <- merge(dt, dt2[, .SD,
                 .SDcols = c("INFO",
                             "FORMAT",
                             "ensembl_transcript_id",
                             "cDNA_change")],
                 by = c("ensembl_transcript_id",
                        "cDNA_change"))

      sdt[, vcf_type := "intersect"]
         return(sdt)

  }

union_vcf <- function(dt, dt2){

    # a function to take the union of annotated variants
      # across VCFs using SnpEff

      rename_vcf_fields(dt, dt2)
      mdt <- merge(dt, dt2[, .SD,
                 .SDcols = c(dt2 %>% names %include% "^INFO",
                             dt2 %>% names %include% "^FORMAT",
                             "ensembl_transcript_id",
                             "cDNA_change")],
                  by = c("ensembl_transcript_id",
                          "cDNA_change"))

      overlaps <- mdt[, paste0(ensembl_transcript_id,
                               cDNA_change) %>%
                        unique]

      sdt <- rbindlist(list(
            mdt,
            dt[!paste0(ensembl_transcript_id, cDNA_change) %chin%
              overlaps],
            dt2[!paste0(ensembl_transcript_id, cDNA_change) %chin%
              overlaps]
              ), use.names = TRUE, fill = TRUE)


      sdt[, vcf_type := "union"]
           return(sdt)

  }

sample_ids <- lapply(ivfdtl, function(dt){
                    dt$sample_id %>% unique
                      }) %>% unique

    # drop vcfs that had no variants before intersect or union attempts to prevent .SDcols bugs
    drop <- lapply(ivfdtl, function(x){

      names(x) %>% length

    }) %>% unlist

    if ((drop == 1) %>% any){

      message(paste(vcfs[which(drop == 1)], "returned no suitable variants and will be excluded from output.\n", sep = " "))

      ivfdtl <- ivfdtl[which(drop != 1)]

      if (length(ivfdtl) == 0) return(NULL)

    }


    # return an intersected data table of variants

      sdt <- parallel::mclapply(sample_ids,

function(sn){

        # find data tables with matching sample names

sdt <- lapply(ivfdtl, function(dt){

             dt[, sample_id %>% .[1]] == sn

            }) %>% unlist

      # merge all data tables with matching sample names

        if (ivfdtl[sdt] %>% length == 1)
          return(ivfdtl[[sdt %>% which]])

        if ((ivfdtl[sdt] %>% length > 1) & intersect == TRUE) {

          x <- ivfdtl[sdt] %>% Reduce(merge_vcf, .)

          if (nrow(x) > 0) return(x)

          if (nrow(x) == 0){

            message(paste("No variants from", sn, "intersect. Returning union."))

            return(ivfdtl[sdt] %>% Reduce(union_vcf, .))

          }

        }


        if ((ivfdtl[sdt] %>% length > 1) & intersect == FALSE)
          return(ivfdtl[sdt] %>% Reduce(union_vcf, .))


      })


  if (!missing(prop_tab)){

    if (length(prop_tab) != length(vcfs)){
      message("Requires one proportions table per vcf file.  Returning vcf data only.")
      sdt %<>% data.table::rbindlist(use.names = TRUE, fill = TRUE)
      return(sdt)
    }

    sdt <- parallel::mclapply(prop_tab %>% seq_along, function(i){

      if (class(prop_tab[i])[1] == "character") pt <- prop_tab[i] %>% rio::import %>% data.table::as.data.table
      if (class(prop_tab[i])[1] == "list") pt <- prop_tab[[i]]
      if (class(pt)[1] == "data.frame") pt %>% data.table::as.data.table

      type <- NULL

      if (!"AF" %chin% names(pt)) type <- "AF"
      if (!"CELLFRACTION" %chin% names(pt)) type <- "CELLFRACTION"

      if (length(type) == 0){

        message("Unable to parse prop_tab, check input table is properly formatted.  Returning with vcf data only.")
        return(sdt)

      }

      if (type == "AF"){

        if (any(!c("#CHROM", "POS", "REF", "ALT", "AF") %chin% names(pt))){
          message("prop_tab is missing columns, see ?garnish_variants.  Returning vcf data only.")
          return(sdt)
        }

        if (dt[, ALT %like% ","] %>% any){
          message("Detected multiple variants per row.  Split allelic fraction input into a single variant per row.  Returning vcf data only.")
          return(sdt)
        }

      # reformat chrom col to account for different reference versions
        pt[, `#CHROM` := `#CHROM` %>% stringr::str_extract("^(chr)?([0-9][0-9]?|[XY]|MT)") %>%
              stringr::str_replace("chr", "")]

        sdt[, `#CHROM` := `#CHROM` %>% stringr::str_extract("^(chr)?([0-9][0-9]?|[XY]|MT)") %>%
                      stringr::str_replace("chr", "")]

        sdt <- merge(sdt, pt, all.x = TRUE, by = c("#CHROM", "POS", "REF", "ALT"))

      }

      if (type == "CELLFRACTION"){

        if (any(!c("chr", "start", "end", "CELLFRACTION") %chin% names(pt))){
          message(paste("prop_tab for", vcfs[i], "is missing columns, see ?garnish_variants.  Returning vcf data only."), sep = " ")
          return(sdt)
        }

        # reformat chrom col to account for different reference versions
        pt[, chr := chr %>% stringr::str_extract("^(chr)?([0-9][0-9]?|[XY]|MT)") %>%
              stringr::str_replace("chr", "")]

        sdt[, `#CHROM` := `#CHROM` %>% stringr::str_extract("^(chr)?([0-9][0-9]?|[XY]|MT)") %>%
                      stringr::str_replace("chr", "")]

        pt %>% data.table::setnames(c("chr", "start"), c("#CHROM", "POS"))

        sdt <- merge(sdt, pt[, .SD %>% unique, .SDcols = c("#CHROM", "POS", "sample_id")],
        all.x = TRUE, by = c("#CHROM", "POS", "sample_id"))

      }

    }) %>% data.table::rbindlist(fill = TRUE, use.names = TRUE)

  }

  if(missing(prop_tab)) sdt %<>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

  if (nrow(sdt) == 0){
    message("All samples returned no suitable variants and will be excluded from output.")
    return(NULL)
  }

# select protein coding variants without NA
sdt %<>%
  .[protein_coding == TRUE &
  !protein_change %>% is.na &
  !effect_type %>% is.na &
  effect_type %like% "insertion|deletion|missense|frameshift"]

  if (nrow(sdt) == 0){
    message("All samples returned no suitable variants and will be excluded from output.")
    return(NULL)
  }

  return(sdt)

}



## ---- garnish_plot
#' Graph `garnish_predictions` results.
#'
#' Plot ADN, CDN, priority, frameshift, and fusion derived `nmers` for class I and class II MHC by `sample_id`.
#'
#' @param input Output from `garnish_predictions` to graph. `input` may be a data table object, list of data tables, or file path to a rio::import-compatible file type. If a list of data tables is provided, unique plots will be generated for each data table.
#'
#' @seealso \code{\link{garnish_predictions}}
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
#' As a side effect: saves graph illustrating the number of neoepitopes in each sample to the working directory. The threshold for inclusion of fusion and frameshift-derived neoepitopes is \eqn{IC_{50}} < 1000nM.
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
      cols <- c("sample_id", names(score_dt) %include% "score")
      score_dt <- score_dt[, .SD, .SDcols = cols]
      score_dt %<>% data.table::melt(id.vars = "sample_id",
                                     measure.vars = names(score_dt) %include% "score")

      score_dt[, MHC := "MHC Class I"] %>%
        .[variable %like% "class_II", MHC := "MHC Class II"] %>%
          .[, variable := variable %>%
            stringr::str_extract("^.*(?=(_class_))")]

      if (score_dt %>% stats::na.omit %>% nrow < 1) return(NULL)

      score_dt <- score_dt[!(MHC == "MHC Class II" & variable == "fitness_scores")]

      score_dt %<>% gplot_names

      if (nrow(score_dt[variable == "classic_top_score"]) != 0){

      g <- ggplot2::ggplot(score_dt[variable == "classic_top_score"],
                           ggplot2::aes(x = sample_id, y = value)) +
           ggplot2::geom_col(ggplot2::aes(fill = variable), col = "black", position = "dodge") +
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
           ggplot2::geom_col(ggplot2::aes(fill = variable), col = "black", position = "dodge") +
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
           ggplot2::geom_col(ggplot2::aes(fill = variable), col = "black", position = "dodge") +
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

  return(NULL)

}


## ---- garnish_targets
#' List the dominant neoepitope sequences for each clone in each sample.
#'
#' Requires that garnish_score was computed from allelic or clonal proportion data.
#'
#' @return A data table with the dominant neoepitope per clone per sample, in rank order of clone frequency.
#'
#' @export garnish_targets
#' @md

garnish_targets <- function(dt){

  if (class(dt)[1] == "character") dt <- dt %>% rio::import %>% data.table::as.data.table

  if (class(dt)[1] == "data.frame") dt %<>% data.table::as.data.table

  dt %<>% data.table::copy

  if (!"clone_id" %chin% names(dt)) stop("No allelic fraction or clonality information in the table.  See ?garnish_variants.")

  dt[, fs := max(fitness_score, na.rm = TRUE), by = c("sample_id", "clone_id")]

  dt <- dt[fitness_score == fs]

  n <- names(dt)[which(names(dt) %chin% c("cDNA_change", "protein_change"))]

  if (length(n) < 1) n <- NULL

  dt <- dt[, .SD %>% unique, .SDcols = c("sample_id", "nmer", "MHC", "external_gene_name", n,
                                        "Consensus_scores", "fitness_score", "iedb_score", "min_DAI", "clone_prop", "clone_id")]

  dt <- dt %>% .[order(sample_id, clone_id)]

  return(dt[Consensus_scores < 1000])

}
