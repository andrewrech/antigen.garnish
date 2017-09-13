library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_variants intersection", {

  # load test data
  "antigen.garnish_example_mutect2.vcf" %>%
      utils::download.file(paste0("http://get.rech.io/", .), .)

  "antigen.garnish_example_strelka.vcf" %>%
    utils::download.file(paste0("http://get.rech.io/", .), .)

  # run test
  dto <- garnish_variants(c("antigen.garnish_example_mutect2.vcf",
                     "antigen.garnish_example_strelka.vcf"))

   testthat::expect_equal(dto %>% length, 21)
   testthat::expect_equal(dto %>% nrow, 4)

    })
