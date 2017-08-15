library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_predictions", {

   if (!check_pred_tools() %>% unlist %>% all) {
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
    }

    # load test data
    dt <- data.table::data.table(
           sample_id = "test",
           ensembl_transcript_id =
           c("ENSMUST00000128119",
             "ENSMUST00000044250",
             "ENSMUST00000018743"),
           cDNA_change = c("c.4988C>T",
                           "c.1114T>G",
                           "c.718T>A"),
           MHC = c("HLA-A*02:01 HLA-DRB1*14:67",
                   "H-2-Kb H-2-IAd",
                   "HLA-A*01:47 HLA-DRB1*03:08")) %>%
      garnish_predictions

        testthat::expect_equal(dt %>% nrow, 552)
        testthat::expect_equal(dt$MHC %>% unique %>% sort,
            c("H-2-IAd",
              "H-2-Kb",
              "HLA-A*01:47",
              "HLA-A*02:01",
              "HLA-DRB1*03:08",
              "HLA-DRB1*14:67"))
          testthat::expect_equal(
                dt[, nmer %>%
                nchar %>%
                unique] %>% sort,
                8:15)

    if (file.exists("antigen.garnish_example.vcf"))
      file.remove("antigen.garnish_example.vcf")
    })
