library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

testthat::test_that("garnish_summary", {

  # load test data
    dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_output.txt")
  # run test
     testthat::expect_equal(garnish_summary(dt),
          data.table::fread("http://get.rech.io/antigen.garnish_example_summary.txt"))

    })