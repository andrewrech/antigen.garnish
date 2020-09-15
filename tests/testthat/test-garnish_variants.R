testthat::test_that("garnish_variants", {
  d <- test_data_dir()

  # load test data
  dt <- "antigen.garnish_example.vcf" %>%
    file.path(d, .) %>%

    # run test
    garnish_variants()

  testthat::expect_equal(dt %>% class() %>% .[1], "data.table")
  testthat::expect_equal(dt$cDNA_change, c("c.4988C>T"))

})
