testthat::test_that("run_mhcnuggets", {

  list.files(pattern = "mhcnuggets.*csv") %>% file.remove
  on.exit(list.files(pattern = "mhcnuggets_.*csv") %>% file.remove)

  skip_pred_tools()
  d <- test_data_dir()

  # load test data
    input <- file.path(d,
      c("antigen.garnish_example_mhcnuggets_input_gru_H-2-KB_412b.csv",
      "antigen.garnish_example_mhcnuggets_input_gru_HLA-A0201_44b3.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_H-2-KB_4632.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_HLA-A0201_4a93.csv"
    ))

    # temporarily copy these to working directory
    file.copy(input, basename(input))

  # run test
    run_mhcnuggets()

   testthat::expect_equal(list.files(pattern = "mhcnuggets_output.*csv") %>% length,
                          4)
   testthat::expect_equal(
                          list.files(pattern = "mhcnuggets_output.*csv") %>%
                            parallel::mclapply(fread) %>%
                            data.table::rbindlist %>%
                            nrow,
                            304
                          )
    })
