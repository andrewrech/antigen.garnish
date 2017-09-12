library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)
library(ShortRead)

testthat::test_that("garnish_jaffa", {

  # load test data
  path <- "antigen.garnish_jaffa_results.csv" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.csv", .)
  fasta_path <- "antigen.garnish_jaffa_results.fasta" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)
  db <- "GRCm38"

  # run test
   dt <- garnish_jaffa(path, db, fasta_path)

   testthat::expect_equal(dt %>% class %>% .[1], "data.table")
   testthat::expect_equal(dt %>% nrow, 471)
   testthat::expect_equal(dt %>% length, 11)
   testthat::expect_equal(dt[, mutant_index %>% sum], 86724)

   if (file.exists("antigen.garnish_jaffa_results.csv"))
   file.remove("antigen.garnish_jaffa_results.csv")
   if (file.exists("antigen.garnish_jaffa_results.fasta"))
     file.remove("antigen.garnish_jaffa_results.fasta")
    })
