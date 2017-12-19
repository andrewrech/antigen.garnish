library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("predictions using all MHC types for all prediction tools", {

    if (parallel::detectCores() < 16) {
    testthat::skip("Skipping long running predictions for all MHC types because ncores < 16")
  }

  # load test data

    MHC <- c("H-2-IAb",
             "H-2-IAd",
             "HLA-DPA1*01:03",
             "HLA-DPA1*02:01",
             "H-2-Db",
             "H-2-Dd",
             "HLA-A*01:01",
             "HLA-A*02:01",
             "HLA-C*14:02",
             "HLA-C*15:02")

      dt <- data.table::data.table(
         MHC = MHC,
         cDNA_change = "c.4988C>T",
         ensembl_transcript_id =
         "ENSMUST00000128119",
         sample_id = "test1")

    # run test
      dto <- dt %>% garnish_predictions

    testthat::expect_true(
        dto[!Consensus_scores %>% is.na, MHC %>% unique %>% length] == 10)

    })
