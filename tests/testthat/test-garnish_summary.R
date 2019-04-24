testthat::test_that("garnish_summary", {

  d <- test_data_dir()

  # load test data
    dt <- data.table::fread(file.path(d, "antigen.garnish_example_affinity_output.txt"))

    dto <- garnish_summary(dt)

  # run test
     testthat::expect_equal(dto$nmers, 154)
     testthat::expect_equal(dto$MHC_binders_class_I, 2)
     testthat::expect_equal(dto$variants, 2)
     testthat::expect_equal(ncol(dto), 19)

     dt <- data.table::fread(file.path(d, "antigen.garnish_PureCN_example_output.txt")) %>%
            .[nmer %chin% c("KLQQLEAAA", "RLICPYMEPI", "LICPYMEPI"), dissimilarity := 0]

    dto2 <- garnish_summary(dt)

    testthat::expect_equal(dto2$nmers, 269)
    testthat::expect_equal(dto2$MHC_binders_class_I, 6)
    testthat::expect_equal(dto2$variants, 4)
    testthat::expect_equal(ncol(dto2), 28)
    testthat::expect_equal(dto2$garnish_score %>% signif(digits = 3), 74900)
    testthat::expect_equal(dto2$fitness_scores_class_I %>% signif(digits = 3), 16.2)

    })
