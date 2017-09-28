library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("run_mhcflurry", {

  if (!check_pred_tools() %>% unlist %>% all){
    testthat::skip("Skipping run_mhcflurry because prediction tools are not in PATH")
    }

  list.files(pattern = "mhcflurry.*csv") %>% file.remove
  on.exit(list.files(pattern = "mhcflurry_.*csv") %>% file.remove)

  # load test data
    input <- c(
      "antigen.garnish_example_mhcflurry_input_ff8dc120-44ee-4e33.csv",
      "antigen.garnish_example_mhcflurry_input_dc7e5ced-6a6e-4f12.csv"
      )

    lapply(input, function(i){
      i %>%
      {utils::download.file(paste0("http://get.rech.io/", .), .)}
      })

    # run test
      run_mhcflurry()

   testthat::expect_equal(list.files(pattern = "mhcflurry_output.*csv") %>% length,
                          2)
   testthat::expect_equal(
                          list.files(pattern = "mhcflurry_output.*csv") %>%
                            parallel::mclapply(fread) %>%
                            data.table::rbindlist %>%
                            nrow,
                            154
                          )
    })