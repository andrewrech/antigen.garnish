testthat::test_that("merge_predictions", {

  list.files(pattern = "antigen.garnish.*_output.*csv") %>% file.remove
  on.exit(list.files(pattern = "antigen.garnish.*_output.*csv") %>% file.remove)

  d <- test_data_dir()

  # load test data
    output <- file.path(d, c(
    "antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5.csv",
    "antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2.csv",
    "antigen.garnish_example_mhcnuggets_output_gru_HLA-A0201_fc5e83cf-8db7-4033.csv",
    "antigen.garnish_example_mhcnuggets_output_lstm_HLA-A0201_2d23d1d7-5622-4930.csv"
  ))

  # put these in home dir for testing
  file.copy(output, basename(output))

  # run test
    dto <- merge_predictions(
      readRDS(file.path(d, "antigen.garnish_netMHC_test_output.RDS")),
      readRDS(file.path(d, "antigen.garnish_merge_predictions_input_dt.RDS"))
    )

  testthat::expect_gt(dto %>% nrow, 250)
  testthat::expect_equal(dto %>% length, 111)

})
