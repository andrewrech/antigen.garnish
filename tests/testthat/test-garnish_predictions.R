library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_predictions vcf", {

   if (!check_pred_tools() %>% unlist %>% all){
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

    })

testthat::test_that("garnish_predictions from Excel file", {

  # load test data
  path <- "antigen.garnish_test_input.xlsx" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_test_input.xlsx", .)


  dt <- garnish_predictions(path = path, predict = FALSE)


  testthat::expect_equal(dt %>% nrow, 552)
  testthat::expect_equal(
    dt[, nmer %>%
         nchar %>%
         unique] %>% sort,
    8:15)

})

testthat::test_that("test predict from jaffa input", {

  if (!check_pred_tools() %>% unlist %>% all){
   testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
   }


  on.exit({
    list.files(pattern = "antigen.garnish_jaffa") %>%
      file.remove})

  # load test data

    path <- "antigen.garnish_jaffa_results.csv" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.csv", .)
    fasta_path <- "antigen.garnish_jaffa_results.fasta" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)

  # get predictions
    dt <- garnish_jaffa(path = path,
                        db = "GRCm38",
                        fasta_path = fasta_path)
    # add mhc values
    dt[, MHC := "H-2-Kb"]

    # predictions
    dt <- garnish_predictions(dt)

    testthat::expect_equal(dt %>% class %>% .[1], "data.table")
    testthat::expect_equal(dt %>% nrow, 36001)
    testthat::expect_equal(dt %>% length, 49)
    testthat::expect_equal(dt[, nmer %>% unique %>% length], 6510)

   })


testthat::test_that("garnish_predictions peptide assemble", {

   if (!check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
    }

  dt <- data.table::data.table(
          sample_id = "test",
          pep_mut = "ATGACTGAATATAAACTTGTGGTA",
          mutant_index = c("7 13 14",  "all", NA),
          MHC = c("H-2-Kb HLA-A*02:01")
                               )

  dto <- garnish_predictions(dt, predict = FALSE)

    testthat::expect_equal(dto %>% length,
                         41)
    testthat::expect_equal(dto %>% nrow,
                         168)
    })

testthat::test_that("garnish_predictions peptide", {

   if (!check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
    }

  dt <- data.table::data.table(
          sample_id = "test",
          pep_mut = "ATGACTGAATATAAACTTGTGGTA",
          mutant_index = "7 13 14",
          MHC = "H-2-Kb"
                               )

  dto <- garnish_predictions(dt)

    testthat::expect_equal(dto %>% length,
                         41)
    testthat::expect_equal(dto %>% nrow,
                         168)
    })

testthat::test_that("garnish_predictions warn on missing pred tools", {

   if (check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping test for warn on missing pred tools because prediction tools are in PATH")
    }

          data.table::data.table(
          sample_id = "test",
          pep_mut = "ATGACTGAATATAAACTTGTGGTA",
          mutant_index = "7 13 14",
          MHC = "H-2-Kb"
                               ) %>%

    {testthat::expect_warning(garnish_predictions(.))}

    })
