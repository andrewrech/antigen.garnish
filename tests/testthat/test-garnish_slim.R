testthat::test_that("garnish_slim", {

  d <- test_data_dir()

  # load test data
  dt <- data.table::fread(file.path(d, "antigen.garnish_full_example_output.txt"))

  # make it slim
  dto <- garnish_slim(dt)

  # run test
  testthat::expect_equal(nrow(dto), 281)

  testthat::expect_equal(max(dto$iedb_score, na.rm = TRUE), 5.94051e-05)

  testthat::expect_equal(nrow(dto[antigen.garnish_class != "Unclassified"]), 3)

  dto <- garnish_slim(dt, slimmer = FALSE)


})
