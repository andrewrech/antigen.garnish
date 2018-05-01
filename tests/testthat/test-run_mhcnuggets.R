testthat::test_that("run_mhcnuggets", {

  list.files(pattern = "mhcnuggets.*csv") %>% file.remove
  on.exit(list.files(pattern = "mhcnuggets_.*csv") %>% file.remove)

  skip_pred_tools()

  # load test data
    input <- file.path(d,
      c("antigen.garnish_example_mhcnuggets_input_gru_H-2-KB_e2bcc9b8-94fa-412b.csv",
      "antigen.garnish_example_mhcnuggets_input_gru_HLA-A0201_9efed71c-199b-44b3.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_H-2-KB_97739391-6066-4632.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_HLA-A0201_e4c2872d-c055-4a93.csv"
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
