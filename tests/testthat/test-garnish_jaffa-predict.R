library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)
library(ShortRead)

testthat::test_that("garnish_jaffa_predict", {

  on.exit({
    list.files(pattern = "antigen.garnish_jaffa") %>%
      file.remove})

  if (!check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
  }

  # load test data

  path <- "antigen.garnish_jaffa_results.csv" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.csv", .)
  fasta_file <- "antigen.garnish_jaffa_results.fasta" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)
  db <- "GRCm38"
  MHCdt <- data.table::data.table(sample_id = c("4662", "Abx7"),
                                  MHC = c("H-2-Kb H-2-Db H-2-IAb", "H-2-Ld H-2-IAd"))

  # run test
  dt <- garnish_jaffa(path, db, MHCdt, fasta_file) %>%
            garnish_predictions

    testthat::expect_equal(dt %>% class %>% .[1], "data.table")
    testthat::expect_equal(dt %>% nrow, 1827040)
    testthat::expect_equal(dt[, nmer %>% unique %>% length], 7776)

   })
