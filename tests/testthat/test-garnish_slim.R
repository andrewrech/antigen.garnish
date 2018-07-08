testthat::test_that("garnish_slim", {

  d <- test_data_dir()

  # load test data
  dt <- data.table::fread(file.path(d, "antigen.garnish_full_example_output.txt"))

  # make it slim
  dto <- garnish_slim(dt)

  # run test
  testthat::expect_equal(nrow(dto), 724)

  unique_dai_scores <- length(unique(dto$DAI))
  testthat::expect_equal(unique_dai_scores, 276)

  max_iedb_score <- max(dto$iedb_score, na.rm = TRUE)
  testthat::expect_equal(max_iedb_score, 5.94051e-05)

  max_NeoRecPotenial <- max(dto$NeoantigenRecognitionPotential, na.rm = TRUE)
  testthat::expect_equal(max_NeoRecPotenial, 2.368751)

  testthat::expect_equal(length(names(dto)), 32)

})
