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

testthat::test_that("garnish_variants", {

  # load test data
    dt <- "antigen.garnish_test_pureCN.vcf" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_test_pureCN.vcf", .) %>%

  # run test
    garnish_variants

  testthat::expect_equal(dt %>% names %>% length, 22)
  testthat::expect_equal(dt[!is.na(CELLFRACTION)] %>% nrow, 29)

    })
    
}


parallel::mclapply(list(vars, m1_CN), test_runner)
