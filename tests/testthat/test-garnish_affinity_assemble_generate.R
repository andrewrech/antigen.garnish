testthat::test_that("garnish_affinity_assemble_generate", {

  list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove
  on.exit(list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove)

  d <- test_data_dir()

  # load test data
    dt <- data.table::fread(file.path(d, "antigen.garnish_example_input.txt")) %>%
     .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                   "H-2-Kb H-2-IAd",
                   "HLA-A*01:47 HLA-DRB1*03:08")]

  # run test
    dto <- garnish_affinity(dt, predict = FALSE, blast = FALSE, fitness = FALSE)

  testthat::expect_equal(dto %>% nrow, 551)
  testthat::expect_equal(dto %>% length, 82)
  testthat::expect_true(dto$cDNA_locs %>% unique %in% c(4988, 718, 1114) %>% all)

})
