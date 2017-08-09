library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("run_netMHC", {

  # load test data
  dt <- "antigen.garnish_example.vcf" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

  # run test
  antigen.garnish::garnish_variants

   testthat::expect_equal(dt %>% class %>% .[1], "data.table")
   testthat::expect_equal(dt$cDNA_change, c("c.4988C>T", "c.1114T>G", "c.718T>A"))

    })
