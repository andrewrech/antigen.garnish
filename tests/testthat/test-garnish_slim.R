testthat::test_that("garnish_slim", {

  d <- test_data_dir()

  # load test data
  dt <- data.table::fread(file.path(d, "antigen.garnish_full_example_output.txt"))

  # make it slim
  dto <- garnish_slim(dt)

  # run test
  testthat::expect_equal(nrow(dto), 281)

  testthat::expect_equal(max(dto$iedb_score, na.rm = TRUE), 5.94051e-05)

  testthat::expect_equal(round(max(dto$NeoantigenRecognitionPotential, na.rm = TRUE), 4), 2.3688)

  testthat::expect_equal(length(names(dto)), 20)

  testthat::expect_equal(nrow(dto[!is.na(antigen.garnish_class)]), 4)

  dto <- garnish_slim(dt, slimmer = FALSE)

  testthat::expect_equal(length(names(dto)), 25)

})
