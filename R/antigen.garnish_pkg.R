#' \pkg{antigen.garnish}: ensemble neoepitope prediction from DNA variants in R.
#'
#' [antigen.garnish](http://neoepitopes.io) is an R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, and gene fusions and performs neoepitope prediction. Input is a VCF file or table of peptides. Output is neoepitopes and a summary of neoepitope load and fitness by sample.
#'
#'Advantages
#'
#'1. **Simplicity**:
#'    - VCF or table input
#'    - summarized neoepitopes for each sample
#'1. **Thoroughness**:
#'    - missense mutations, insertions, deletions, and gene fusions
#'    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i).
#'    - filter against all known normal proteins/immunogenic epitopes and then rank by fitness
#'1. **Speed**:
#'    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
#'    - on an Amazon Web Services `m4.16xlarge` EC2 instance, 20,000 consensus predictions using 100+ MHC types in under 5 minutes
#' @section Manifest:
#'* `garnish_variants`: process missense / indel VCF variants from [SnpEff](http://snpeff.sourceforge.net/)
#'* `garnish_jaffa`: process gene fusions from [JAFFA](https://github.com/Oshlack/JAFFA)
#'* `garnish_predictions`: perform ensemble neoepitope prediction
#'* `garnish_summary`: summarize and rank results
#'* `garnish_plot`: generate summary plots
#'* `list_mhc`: list all supported MHC allele syntax
#' @section Bug, issues, feedback:
#' Please report bugs and issues and provide feedback via [Github](https://github.com/andrewrech/antigen.garnish/issues) or to [andrewrech\@gmail.com](mailto:andrewrech\@gmail.com).
#' @docType package
#' @name antigen.garnish
#' @name BiocInstaller
#' @import colorspace
#' @import data.table
#' @import dt.inflix
#' @import ggplot2
#' @import parallel
#' @import stringr
#' @import testthat
#' @importFrom rio import
#' @importFrom ShortRead sread readFasta
#' @importFrom Biostrings DNAString translate
#' @importFrom biomaRt useMart getBM getSequence
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom Rdpack reprompt
#' @importFrom stats na.omit t.test
#' @importFrom stringi stri_detect_fixed
#' @importFrom tidyr separate_rows
#' @importFrom utils download.file
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#' @md
NULL
