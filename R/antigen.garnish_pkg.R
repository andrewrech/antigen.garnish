#' \pkg{antigen.garnish}: ensemble neoepitope prediction from DNA variants in R.
#'
#' An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, and gene fusions and performs ensemble neoepitope prediction using 7 algorithms. Input is a VCF file or table of peptides. Output is ranked neoepitopes and a summary of neoepitope load and fitness by sample. Neoepitopes are ranked by MHC I/II binding affinity, clonality, RNA expression, dissimilarity to the normal peptidome, and similarity to known immunogenic antigens. [More information.](http://antigen-garnish-presentation.s3-website-us-east-1.amazonaws.com)
#'
#'Advantages
#'
#'1. **Thoroughness**:
#'    - missense mutations, insertions, deletions, and gene fusions
#'    - human and mouse
#'    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
#'    - ranked by
#'    	+ MHC I/II binding affinity
#'    	+ clonality
#'    	+ RNA expression
#'    	+ dissimilarity to the normal peptidome
#'    	+ similarity to known immunogenic antigens
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
#'* `garnish.antigens`: list top neoepitope sequences, sources, and properties for each [clone](https://github.com/lima1/PureCN) and sample
#'* `garnish_plot`: generate summary plots
#' @section Bug, issues, feedback:
#' Please report bugs and issues and provide feedback via [Github](https://github.com/andrewrech/antigen.garnish/issues) or to [andrewrech\@gmail.com](mailto:andrewrech\@gmail.com).
#' @docType package
#' @name antigen.garnish
#' @name BiocInstaller
#' @import colorspace
#' @import data.table
#' @import dt.inflix
#' @import ggplot2
#' @import knitr
#' @import mclust
#' @import mclust
#' @import parallel
#' @import stringr
#' @import testthat
#' @importFrom Biostrings DNAString translate
#' @importFrom Rdpack reprompt
#' @importFrom ShortRead sread readFasta
#' @importFrom biomaRt useMart getBM getSequence
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom rio import
#' @importFrom stats na.omit t.test
#' @importFrom stringi stri_detect_fixed
#' @importFrom tidyr separate_rows
#' @importFrom utils download.file packageVersion
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#' @md
NULL