testthat::test_that("garnish_affinity_assemble_generate", {

  list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove
  on.exit(list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove)

  d <- test_data_dir()

  # load test data
    dt <- data.table::fread(file.path(d, "antigen.garnish_example_input.txt")) %>%
     .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67")]

  # run test
    dto <- garnish_affinity(dt, predict = FALSE, blast = FALSE, fitness = FALSE)

  testthat::expect_equal(dto %>% nrow, 183)
  testthat::expect_equal(dto %>% length, 78)
  testthat::expect_true(identical(dto$cDNA_locs %>% unique, as.integer(4988)))

})
