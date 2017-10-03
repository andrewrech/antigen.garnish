library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)
library(ggplot2)
library(cowplot)

testthat::test_that("garnish_plot", {
  
  on.exit(list.files(pattern = "(antigen.garnish.*\\.pdf)|(antigen.garnish_example_jaffa)") %>% file.remove)
  
  # load test data
  dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_jaffa_output.txt")
  
  jdt <- dt[pep_type == "fus"] %>% .[, fus_tx := "AAAGGGCCCCTTT"]
  
  dt <- dt[pep_type != "fus"]
  
  g <- garnish_plot(dt = dt, jdt = jdt)
  
  # run test
  testthat::expect_equal(list.files(pattern = "\\.pdf") %>% length, 3)
  testthat::expect_equal(length(g), 3)
  testthat::expect_equal(class(g[[1]])[1], "gg")
  
})