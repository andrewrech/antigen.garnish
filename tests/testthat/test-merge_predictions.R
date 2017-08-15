library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_predictions", {

# load test data and run test
  # requires mhcflurry output in disk
    "antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5-94e3-c8fcfe407669.csv" %>%
     utils::download.file(
          paste0("http://get.rech.io/", .), .)

    "antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2-88f9-0628be05c1ea.csv" %>%
     utils::download.file(
          paste0("http://get.rech.io/", .), .)

    dto <- merge_predictions(
            readRDS(gzcon(url("http://get.rech.io/antigen.garnish_netMHC_test_output.RDS"))),
            readRDS(gzcon(url("http://get.rech.io/antigen.garnish_merge_predictions_input_dt.RDS")))
                   )

    testthat::expect_equal(dto %>% nrow, 552)
    testthat::expect_equal(dto %>% length, 107)

    if (file.exists("antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5-94e3-c8fcfe407669.csv"))
    file.remove("antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5-94e3-c8fcfe407669.csv")


    if (file.exists("antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2-88f9-0628be05c1ea.csv"))
    file.remove("antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2-88f9-0628be05c1ea.csv")

    })
