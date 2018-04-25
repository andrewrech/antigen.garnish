vars <- function(){

testthat::test_that("garnish_variants", {

  # load test data
    dt <- "antigen.garnish_example.vcf" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

  # run test
    garnish_variants

  testthat::expect_equal(dt %>% class %>% .[1], "data.table")
  testthat::expect_equal(dt$cDNA_change, c("c.4988C>T", "c.1114T>G", "c.718T>A"))

    })

}

m1_CN <- function(){

testthat::test_that("garnish_variants with pureCN", {

  # load test data
    dt <- "antigen.garnish_test_pureCN.vcf" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_test_pureCN.vcf", .) %>%

  # run test
    garnish_variants

  testthat::expect_equal(dt %>% names %>% length, 22)
  testthat::expect_equal(dt[!is.na(CELLFRACTION)] %>% nrow, 29)

    })

}

m1_AF <- function(){

testthat::test_that("garnish_variants with AF/prop_tab", {

  # load test data
  dt <- "antigen.garnish_example.vcf" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .)


    # make a mock table
    pt <- vcfR::read.vcfR(dt) %>% .@fix %>% data.table::as.data.table

    pt <- pt[, .SD %>% unique, .SDcols = c("CHROM", "POS", "REF", "ALT")] %>%
            .[, AF := 0.5] %T>%
            data.table::fwrite("antigen.garnish_prop_AF.txt", sep = "\t", quote = FALSE, row.names = FALSE)

    pt[, .SD %>% unique, .SDcols = c("CHROM", "POS")] %>%
    .[, CELLFRACTION := c(0.25, 0.3, .4)] %>%
    .[, end := POS] %>%
    data.table::setnames(c("CHROM", "POS"), c("chr", "start")) %>%
    data.table::fwrite("antigen.garnish_prop_CF.txt", sep = "\t", quote = FALSE, row.names = FALSE)

  # run test
    dt <- garnish_variants(dt, prop_tab = "antigen.garnish_prop_AF.txt")

	  testthat::expect_equal(dt %>% names %>% length, 22)
	  testthat::expect_equal(dt[is.na(AF)] %>% nrow, 0)

	  dt <- garnish_variants(dt, prop_tab = "antigen.garnish_prop_CF.txt")

	  testthat::expect_equal(dt %>% names %>% length, 22)
	  testthat::expect_equal(dt[, CELLFRACTION], c(0.3, 0.4, 0.25))

    })

}


parallel::mclapply(list(vars, m1_CN, m1_AF), test_runner)
