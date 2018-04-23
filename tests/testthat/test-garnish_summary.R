testthat::test_that("garnish_summary", {

  # load test data
    dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_output.txt")

    dto <- garnish_summary(dt)

  # run test
     testthat::expect_equal(dto$nmers, 154)
     testthat::expect_equal(dto$mhc_binders_class_I, 3)
     testthat::expect_equal(dto$variants, 2)
     testthat::expect_equal(ncol(dto), 21)

     dt <- data.table::fread("http://get.rech.io/antigen.garnish_pureCN_example_output.txt")

    dto2 <- garnish_summary(dt)

    testthat::expect_equal(dto2$nmers, 269)
    testthat::expect_equal(dto2$mhc_binders_class_I, 8)
    testthat::expect_equal(dto2$variants, 4)
    testthat::expect_equal(ncol(dto2), 24)
    testthat::expect_equal(dto2$garnish_score %>% signif(digits = 3), 74900)
    testthat::expect_equal(dto2$fitness_scores_class_I %>% signif(digits = 3), 16.2)

    })
