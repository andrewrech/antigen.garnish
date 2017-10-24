library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_predictions README example from VCF", {

   if (!check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
    }

    # load test data
      dt <- "antigen.garnish_example.vcf" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

    # run test
      garnish_variants %>%
        .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                     "H-2-Kb H-2-IAd",
                    "HLA-A*01:47 HLA-DRB1*03:08")] %>%
      garnish_predictions %>%
      garnish_summary

    testthat::expect_equivalent(
          dt,
          data.table::data.table(sample_id = "normal_tumor.bam",
                         priority_neos_class_I = 0L,
                         priority_neos_class_II = 0L,
                         classic_neos_class_I = 0L,
                         classic_neos_class_II = 0L,
                         classic_top_score_class_I = 0.0106647742091908,
                         classic_top_score_class_II = 0.0010342955988725,
                         alt_neos_class_I = 3L,
                         alt_neos_class_II = 0L,
                         alt_top_score_class_I = 33.424955713995,
                         alt_top_score_class_II = 3.00997014163189,
                         fs_neos_class_I = 0L,
                         fs_neos_class_II = 0L,
                         fusion_neos_class_I = 0L,
                         fusion_neos_class_II = 0L,
                         mhc_binders_class_I = 7L,
                         mhc_binders_class_II = 6L,
                         predictions = 275L,
                         nmers = 275L,
                         variants = 3L,
                         transcripts = 3L)
          )
    })

testthat::test_that("garnish_predictions from transcripts", {

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
      # run test
        garnish_predictions

      testthat::expect_equal(dt %>% nrow, 551)
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

  # run test
  dt <- garnish_predictions(path = path, predict = FALSE)


  testthat::expect_equal(dt %>% nrow, 551)
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

  # load test data
    dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_input_jaffa.csv")

    dt[, MHC := "H-2-Kb"]

  # run test
    dt <- garnish_predictions(dt)

  testthat::expect_equal(dt %>% class %>% .[1], "data.table")
  testthat::expect_equal(dt %>% nrow, 1071)
  testthat::expect_equal(dt[, nmer %>% unique %>% length], 140)

   })

testthat::test_that("garnish_predictions assemble from peptides", {

  # load test data
    dt <- data.table::data.table(
            sample_id = "test",
            pep_mut = "AAVMILKWTFRAKINGDEHKRDEF",
            mutant_index = c("7 13 14",  "all", NA),
            MHC = c("H-2-Kb HLA-A*02:01")
                                 )
  # run test data
    dto <- garnish_predictions(dt, predict = FALSE)

  testthat::expect_equal(dto %>% length,
                         13)
  testthat::expect_equal(dto %>% nrow,
                         880)
    })

testthat::test_that("garnish_predictions peptide", {

   if (!check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
    }

  # load test data
    dt <- data.table::data.table(
            sample_id = "test",
            pep_mut = "AAVMILKWTFRAKINGDEHKRDEF",
            mutant_index = "7 13 14",
            MHC = "H-2-Kb"
                                 )
  # run test
    dto <- garnish_predictions(dt)

  testthat::expect_equal(dto %>% length,
                         42)
  testthat::expect_equal(dto %>% nrow,
                         168)
    })

testthat::test_that("garnish_predictions warning on missing prediction tools", {

   if (check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping test for warning on missing pred tools because prediction tools are in PATH")
    }
   # load test data
      data.table::data.table(
      sample_id = "test",
      pep_mut = "ATGACTGAATATAAACTTGTGGTA",
      mutant_index = "7 13 14",
      MHC = "H-2-Kb"
                           ) %>%
   # run test
    {
      testthat::expect_warning(garnish_predictions(.))
    }

    })
