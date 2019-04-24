testthat::test_that("merge_predictions", {

  list.files(pattern = "antigen.garnish.*_output.*csv") %>% file.remove
  on.exit(list.files(pattern = "antigen.garnish.*_output.*csv") %>% file.remove)

  d <- test_data_dir()

  # load test data
    output <- file.path(d, c(
    "a.g_ex_mhcflurry_output_42c5.csv",
    "a.g_ex_mhcflurry_output_47f2.csv",
    "a.g_ex_mhcnuggets_output_gru_HLA-A0201_4033.csv",
    "a.g_ex_mhcnuggets_output_lstm_HLA-A0201_4930.csv"
  ))

  # put these in home dir for testing
  file.copy(output, basename(output))

  # run test

  dto <- merge_predictions(
   file.path(d, "antigen.garnish_netMHC_test_output.RDS") %>%
      readRDS,
  file.path(d, "antigen.garnish_merge_predictions_input_dt.RDS") %>%
      readRDS
 )

  testthat::expect_gt(dto %>% nrow, 250)
  testthat::expect_equal(dto %>% length, 111)

})
