testthat::test_that("garnish_affinity_assemble_generate", {

  list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove
  on.exit(list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove)

  d <- test_data_dir()

  # load test data
    dt <- data.table::fread(file.path(d, "antigen.garnish_example_input.txt")) %>%
     .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67")]

  # run test
    dto <- garnish_affinity(dt, predict = FALSE, blast = FALSE)

  testthat::expect_equal(dto %>% nrow, 183)
  testthat::expect_equal(dto %>% length, 78)
  testthat::expect_true(identical(dto$cDNA_locs %>% unique, as.integer(4988)))

  # test edge cases
  dt <- data.table::fread(file.path(d, "antigen.garnish_example_frameshifts.txt"))

  dto <- garnish_affinity(dt, generate = FALSE, predict = FALSE,
    blast = FALSE, remove_wt = FALSE)

  testthat::expect_equal(dto[, mutant_index],
    c("2790", "2789", "2790", "8", "9 10", "10",
      "10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28",
    "10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28"))

})
