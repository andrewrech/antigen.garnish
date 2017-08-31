library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("list_mhc", {

  testthat::expect_equal(
    list_mhc() %>% length,
    3874)

    })
