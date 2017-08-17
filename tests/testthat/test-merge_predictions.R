library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_predictions", {
  
  # load test data and run test
  # requires mhcflurry output in disk
  file.copy(from = system.file("extdata",
                               "ag_ex_mhcnuggets_output_gru_HLA-A0201_fc5e83cf-8db7-4033.csv",
                               package = "antigen.garnish"), to = getwd())
  
  file.copy(from = system.file("extdata",
                               "ag_ex_mhcnuggets_output_lstm_HLA-A0201_2d23d1d7-5622-4930.csv",
                               package = "antigen.garnish"), to = getwd())
  
  
  file.copy(from = system.file("extdata",
                               "antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2.csv",
                               package = "antigen.garnish"), to = getwd())
  
  file.copy(from = system.file("extdata",
                               "antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5.csv",
                               package = "antigen.garnish"), to = getwd())
  
  dto <- merge_predictions(
    readRDS(gzcon(url("http://get.rech.io/antigen.garnish_netMHC_test_output.RDS"))),
    readRDS(gzcon(url("http://get.rech.io/antigen.garnish_merge_predictions_input_dt.RDS")))
  )
  
  testthat::expect_equal(dto %>% nrow, 552)
  testthat::expect_equal(dto %>% length, 110)
  
  if (file.exists("ag_ex_mhcnuggets_output_gru_HLA-A0201_fc5e83cf-8db7-4033.csv"))
    file.remove("ag_ex_mhcnuggets_output_gru_HLA-A0201_fc5e83cf-8db7-4033.csv")
  
  
  if (file.exists("ag_ex_mhcnuggets_output_lstm_HLA-A0201_2d23d1d7-5622-4930.csv"))
    file.remove("ag_ex_mhcnuggets_output_lstm_HLA-A0201_2d23d1d7-5622-4930.csv")
  
  if (file.exists("antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2.csv"))
    file.remove("antigen.garnish_example_mhcflurry_output_68e72a3d-3504-47f2.csv")
  
  if (file.exists("antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5.csv"))
    file.remove("antigen.garnish_example_mhcflurry_output_6e7d571f-8b47-42c5.csv")
  
})

