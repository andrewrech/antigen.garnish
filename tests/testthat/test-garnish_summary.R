library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_summary", {

  # load test data
    dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_output.txt")

    dto <- garnish_summary(dt)

  # run test
     testthat::expect_equal(dto$nmers, 154)
     testthat::expect_equal(dto$mhc_binders, 7)
     testthat::expect_equal(dto$variants, 2)
     testthat::expect_equal(dto %>% names,
      c("sample_id",
        "priority_neos",
        "classic_neos",
        "classic_top_score",
        "alt_neos",
        "alt_top_score",
        "mhc_binders",
        "variants",
        "transcripts",
        "predictions",
        "nmers"))
    })