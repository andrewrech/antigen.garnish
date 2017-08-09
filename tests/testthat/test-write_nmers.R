library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

test_that("write_nmers", {

  # load test data
  dt <-

  write_nmers

   testthat::expect_equal(dt %>% class %>% .[1], "data.table")
   testthat::expect_equal(dt$cDNA_change, c("c.4988C>T", "c.1114T>G", "c.718T>A"))

    })
