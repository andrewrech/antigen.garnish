library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("cat_tables", {
  
  on.exit(list.files(pattern = "antigen.garnish_example_jaffa_output") %>% file.remove)
  
  # load test data
  dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_jaffa_output.txt")
  
  jdt <- dt[pep_type == "fus"]
  
  dt <- dt[pep_type != "fus"]
  
  dt <- cat_tables(jdt, dt)
  
  # run test
  testthat::expect_equal(dt %>% nrow, 100)
  testthat::expect_equal(dt[, nmer %>% unique %>% length], 95)
  testthat::expect_equal(dt[pep_type != "fus"] %>% nrow, 12)
  
})
