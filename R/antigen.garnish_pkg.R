#' \pkg{antigen.garnish}: ensemble neoantigen prediction from DNA variants in R.
#'
#'An R package for [neoantigen](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, or RNASeq-derived gene fusions and performs ensemble neoantigen prediction using 7 algorithms. Input is a VCF file, [JAFFA](https://github.com/Oshlack/JAFFA) output, or table of peptides or transcripts. Outputs are ranked and summarized by sample. Neoantigens are ranked by MHC I/II binding affinity, clonality, RNA expression, similarity to known immunogenic antigens, and dissimilarity to the normal peptidome.
#'
#'Advantages
#'
#'1. **Thoroughness**:
#'  * missense mutations, insertions, deletions, and gene fusions
#'  * human and mouse
#'  * ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
#'  * ranked by
#'    * MHC I/II binding affinity
#'    * clonality
#'    * RNA expression
#'    * similarity to known immunogenic antigens
#'    * dissimilarity to the normal peptidome
#'2. **Speed and simplicity**:
#'  * 1000 variants are ranked in a single step in less than five minutes
#'  * parallelized using [`parallel::mclapply`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html) and [data.table::setDTthreads](https://github.com/Rdatatable/data.table/wiki), see respective links for information on setting multicore usage
#'3. **Integration with R/Bioconductor**
#'  * upstream/VCF processing
#'  * exploratory data analysis, visualization
#'
#'[Package documentation](https://neoantigens.rech.io/reference/index.html) ([pdf](https://get.rech.io/antigen.garnish.pdf))
#'
#'We welcome contributions and feedback via [Github](https://github.com/immune-health/antigen.garnish/issues) or [email](mailto:leepr@upenn.edu).
#'
#' @seealso \code{\link{garnish_variants}}
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_summary}}
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
#' @importFrom testthat expect_equal expect_equivalent expect_gt expect_true skip succeed test_that
#' @importFrom Biostrings DNAString translate pairwiseAlignment
#' @importFrom Rdpack reprompt
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom rio import
#' @importFrom stats na.omit t.test ecdf quantile
#' @importFrom stringi stri_detect_fixed
#' @importFrom tidyr separate_rows
#' @importFrom utils download.file packageVersion installed.packages
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#'
#' @md
NULL
