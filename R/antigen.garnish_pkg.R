#' \pkg{antigen.garnish}: ensemble neoepitope prediction from DNA variants in R.
#'
#' [antigen.garnish](http://neoepitopes.io) is an R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA point mutations, insertions, and deletions and performs neoepitope prediction. Input is a VCF file or table of peptides. Output is neoepitopes and a summary by sample.
#'
#'Advantages
#'
#'1. **Simplicity**:
#'    - VCF or table input
#'    - summarized neoepitopes for each sample
#'1. **Thoroughness**:
#'    - missense mutations and frameshifts
#'    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i).
#'1. **Speed**:
#'    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
#'    - on an Amazon Web Services `m4.16xlarge` EC2 instance, 20,000 consensus predictions using 100+ MHC types in under 5 minutes
#' @section Manifest:
#' * `garnish_variants`: process variants from [SnpEff](http://snpeff.sourceforge.net/)
#' * `garnish_predictions`: perform ensemble neoepitope prediction
#' * `garnish_summary`: summarize neoepitope prediction
#' @docType package
#' @name antigen.garnish
#' @import colorspace
#' @import data.table
#' @import dt.inflix
#' @import parallel
#' @import stringr
#' @import testthat
#' @importFrom rio import
#' @importFrom Biostrings DNAString translate
#' @importFrom biomaRt useMart getBM getSequence
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom stats na.omit t.test
#' @importFrom stringi stri_detect_fixed
#' @importFrom tidyr separate_rows
#' @importFrom utils download.file
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#' @md
NULL
