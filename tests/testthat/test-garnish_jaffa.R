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
  fasta_file <- "antigen.garnish_jaffa_results.fasta" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)
  db <- "GRCm38"
  MHCdt <- data.table::data.table(sample_id = c("4662", "Abx7"),
                                  MHC = c("H-2-Kb H-2-Db H-2-IAb", "H-2-Ld H-2-IAd"))

  # run test
<<<<<<< HEAD
   dt <- garnish_jaffa(path, db, MHCdt, fasta_file)
=======
  dt <- garnish_jaffa(path, db, MHCdt, fasta_file)
>>>>>>> 88daccfe250a85b69b6441dc33c86782fda8f44b

   testthat::expect_equal(dt %>% class %>% .[1], "data.table")
   testthat::expect_equal(dt %>% nrow, 471)
   testthat::expect_equal(dt %>% length, 12)
   testthat::expect_equal(dt[, mutant_index %>% sum], 86724)

   if (file.exists("antigen.garnish_jaffa_results.csv"))
   file.remove("antigen.garnish_jaffa_results.csv")
   if (file.exists("antigen.garnish_jaffa_results.fasta"))
     file.remove("antigen.garnish_jaffa_results.fasta")
    })
