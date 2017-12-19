library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_jaffa", {

    list.files(pattern = "antigen.garnish_jaffa") %>% file.remove
    on.exit(list.files(pattern = "antigen.garnish_jaffa") %>% file.remove)

  # load test data
    path <- "antigen.garnish_jaffa_results.csv" %T>%
        utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.csv", .)
    fasta_path <- "antigen.garnish_jaffa_results.fasta" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)

  # run test
   dt <- garnish_jaffa(path = path, db = "GRCm38", fasta_path = fasta_path)

  testthat::expect_equal(dt %>% nrow, 15)
  testthat::expect_equal(dt %>% length, 11)
  testthat::expect_equal(dt[, mutant_index %>% sum], 662)
    })
