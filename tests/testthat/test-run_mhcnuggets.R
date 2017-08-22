library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("run_netMHC", {

   if (!check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
  }


    # load test data
    input <- c(
      "antigen.garnish_example_mhcnuggets_input_gru_H-2-KB_dbd00aa2-95de-426d.csv",
      "antigen.garnish_example_mhcnuggets_input_gru_HLA-A0301_aae53bd2-3f2e-4198.csv",
      "antigen.garnish_example_mhcnuggets_input_gru_HLA-B3901_d44c0596-ba5d-4904.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_HLA-A1101_4a35158a-1fc2-42a1.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_HLA-A3201_28a9297e-58e3-48b6.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_HLA-B1801_42078519-6467-42e1.csv"
      )

    on.exit({
      suppressWarnings(file.remove(input, showWarnings = FALSE) %>% invisible)
      })

    parallel::mclapply(input, function(i){
      i %>%
      {utils::download.file(paste0("http://get.rech.io/", .), .)}
      })

    # run test
      run_mhcnuggets()

      dt <- list.files(pattern = "antigen.garnish_example_mhcnuggets_output.*") %>%
        lapply(fread) %>% rbindlist

   testthat::expect_equal(dt %>% nrow, 224)
   testthat::expect_equal(dt %>% length, 2)
    })