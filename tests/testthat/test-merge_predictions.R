testthat::test_that("antigen.garnish:::merge_predictions", {
  list.files(pattern = "antigen.garnish.*_output.*csv") %>% file.remove()
  on.exit(list.files(pattern = "antigen.garnish.*_output.*csv") %>% file.remove())

  d <- test_data_dir()

  # load test data
  output <- file.path(d, c(
    "a.g_ex_mhcflurry_output_42c5.csv",
    "a.g_ex_mhcflurry_output_47f2.csv"
  ))

  # put these in home dir for testing
  file.copy(output, basename(output))

  # run test

  dto <- antigen.garnish:::merge_predictions(
    file.path(d, "antigen.garnish_netMHC_test_output.RDS") %>%
      readRDS(),
    file.path(d, "antigen.garnish_merge_predictions_input_dt.RDS") %>%
      readRDS()
  )

  testthat::expect_gt(dto %>% nrow(), 250)
  testthat::expect_equal(dto %>% length(), 106)
})
