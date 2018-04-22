#' \pkg{antigen.garnish}: ensemble neoepitope prediction from DNA variants in R.
#'
#' [antigen.garnish](http://neoepitopes.io) is an R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, and gene fusions and performs neoepitope prediction. Input is a VCF file or table of peptides. Output is neoepitopes and a summary of neoepitope load and fitness by sample.
#'
#'Advantages
#'
#'1. **Thoroughness**:
#'    - missense mutations, insertions, deletions, and gene fusions
#'    - human and mouse
#'    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
#'    - rank by dissimilarity to the normal peptidome / similarity to known immunogenic antigens
#'1. **Speed and simplicity**:
#'    - 1000 variants are ranked in a single step in less than five minutes
#'1. **Integration with R/Bioconductor**
#'    - upstream/VCF processing
#'    - exploratory data analysis, visualization
#' @section Manifest:
#'* `garnish_variants`: process missense / indel VCF variants from [SnpEff](http://snpeff.sourceforge.net/)
#'* `garnish_jaffa`: process gene fusions from [JAFFA](https://github.com/Oshlack/JAFFA)
#'* `garnish_predictions`: perform ensemble neoepitope prediction
#'* `garnish_summary`: summarize and rank results
#'* `garnish_plot`: generate summary plots
#'* `garnish_targets`: list dominant neoepitope sequences, sources, and properties per clone per sample
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
#' @import mclust
#' @import parallel
#' @import stringr
#' @import testthat
#' @import knitr
#' @importFrom rio import
#' @importFrom ShortRead sread readFasta
#' @importFrom Biostrings DNAString translate
#' @importFrom biomaRt useMart getBM getSequence
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom Rdpack reprompt
#' @importFrom stats na.omit t.test
#' @importFrom stringi stri_detect_fixed
#' @importFrom tidyr separate_rows
#' @importFrom utils download.file packageVersion
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#' @md
NULL