testthat::test_that("garnish_summary", {

  # load test data
    dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_output.txt")

    dto <- garnish_summary(dt)

  # run test
     testthat::expect_equal(dto$nmers, 154)
     testthat::expect_equal(dto$mhc_binders_class_I, 3)
     testthat::expect_equal(dto$variants, 2)

     dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_output.txt") %>%
      garnish_fitness

    dto <- garnish_summary(dt)

    testthat::expect_equal(dto$nmers, 154)
    testthat::expect_equal(dto$mhc_binders_class_I, 7)
    testthat::expect_equal(dto$variants, 2)

    })
