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
     testthat::expect_equal(dto$mhc_binders_class_I, 7)
     testthat::expect_equal(dto$variants, 2)
     testthat::expect_equal(dto %>% names,
      c("sample_id",
        "priority_neos_class_I",
        "priority_neos_class_II",
        "classic_neos_class_I",
        "classic_neos_class_II",
        "classic_top_score_class_I",
        "classic_top_score_class_II",
        "alt_neos_class_I",
        "alt_neos_class_II",
        "alt_top_score_class_I",
        "alt_top_score_class_II",
        "fs_neos_class_I",
        "fs_neos_class_II",
        "fusion_neos_class_I",
        "fusion_neos_class_II",
        "mhc_binders_class_I",
        "mhc_binders_class_II",
        "predictions",
        "nmers",
        "variants",
        "transcripts")
)
    })