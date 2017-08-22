library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("merge_predictions", {

  # load test data
  # requires mhcflurry and mhcnuggets output on disk

    output <- c(
    "antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5.csv",
    "antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2.csv",
    "antigen.garnish_example_mhcnuggets_output_gru_HLA-A0201_fc5e83cf-8db7-4033.csv",
    "antigen.garnish_example_mhcnuggets_output_lstm_HLA-A0201_2d23d1d7-5622-4930.csv"
                )

    on.exit({
      suppressWarnings(file.remove(output, showWarnings = FALSE) %>% invisible)
    })

parallel::mclapply(output, function(i){
        i %>%
            {utils::download.file(paste0("http://get.rech.io/", .), .)}
      })

  # run test
    dto <- merge_predictions(
      readRDS(gzcon(url("http://get.rech.io/antigen.garnish_netMHC_test_output.RDS"))),
      readRDS(gzcon(url("http://get.rech.io/antigen.garnish_merge_predictions_input_dt.RDS")))
    )

  testthat::expect_equal(dto %>% nrow, 552)
  testthat::expect_equal(dto %>% length, 111)

})


saveRDS(dt, "antigen.garnish_merge_predictions_input_dt.RDS")