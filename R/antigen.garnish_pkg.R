#' \pkg{antigen.garnish}: ensemble neoantigen prediction from DNA variants in R.
#'
#' Github: https://github.com/immune-health/antigen.garnish
#' Documentation: https://neoantigens.rech.io/
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
#' @importFrom stats na.omit t.test ecdf quantile
#' @importFrom tidyr separate_rows
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#'
#' @md
NULL
