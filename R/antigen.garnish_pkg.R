#' \pkg{antigen.garnish}: ensemble neoantigen prediction from DNA variants in R.
#'
#' Ensemble tumor neoantigen prediction from complex variants. Immunogenicity filtering based on the Tumor Neoantigen Selection Alliance (TESLA).
#'
#' Advantages
#'
#' 1. **Thoroughness**:
#'   - missense mutations, insertions, or deletions
#'   - human and mouse
#'   - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets-2.0), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
#'   - ranked by
#'     - MHC I/II binding affinity
#'     - clonality
#'     - RNA expression
#'     - similarity to known immunogenic antigens
#'     - dissimilarity to the normal peptidome
#' 2. **Speed and simplicity**:
#'   - 1000 variants are ranked in a single step in less than five minutes
#'   - parallelized using [`parallel::mclapply`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html), [`data.table::setDTthreads`](https://github.com/Rdatatable/data.table/wiki), and [GNU parallel](https://www.gnu.org/software/parallel/)
#' 3. **Integration with R/Bioconductor**
#'   - upstream/VCF processing
#'   - exploratory data analysis, visualization
#'
#' [Package documentation](https://neoantigens.rech.io/reference/index.html) ([pdf](https://get.rech.io/antigen.garnish.pdf))
#'
#' We welcome contributions and feedback via [Github](https://github.com/immune-health/antigen.garnish/issues) or [email](mailto:leepr@upenn.edu).
#'
#' @seealso \code{\link{garnish_variants}}
#' @seealso \code{\link{garnish_affinity}}
#' @seealso \code{\link{garnish_antigens}}
#' @docType package
#' @name antigen.garnish
#' @import data.table
#' @import mclust
#' @import parallel
#' @import stringr
#' @importFrom testthat expect_equal expect_equivalent expect_gt expect_true skip succeed test_that
#' @importFrom Biostrings AAStringSet AMINO_ACID_CODE writeXStringSet DNAString translate pairwiseAlignment
#' @importFrom Rdpack reprompt
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom purrr partial
#' @importFrom rio import
#' @importFrom stats na.omit t.test ecdf quantile
#' @importFrom tidyr separate_rows
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#'
#' @md
NULL
