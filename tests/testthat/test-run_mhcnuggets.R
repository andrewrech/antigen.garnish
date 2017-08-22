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
      "antigen.garnish_example_mhcnuggets_input_gru_H-2-KB_e2bcc9b8-94fa-412b.csv",
      "antigen.garnish_example_mhcnuggets_input_gru_HLA-A0201_9efed71c-199b-44b3.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_H-2-KB_97739391-6066-4632.csv",
      "antigen.garnish_example_mhcnuggets_input_lstm_HLA-A0201_e4c2872d-c055-4a93.csv"
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

   testthat::expect_equal(dt %>% nrow, 304)
   testthat::expect_equal(dt %>% length, 2)
    })