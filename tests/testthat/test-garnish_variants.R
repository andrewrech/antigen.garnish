library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

test_that("garnish_variants", {
browser()
dt <-
    # load an example VCF
    system.file("extdata",
          "antigen.garnish_example.vcf",
          package = "antigen.garnish") %>%

    # extract variants
    antigen.garnish::garnish_variants(.)

testthat::expect_equal(class(dt[[1]])[1], "data.table")
testthat::expect_equal(dt[[3]]$aa_mutation,
                       c("W240R", "S372A", "A1663V"))
    })
