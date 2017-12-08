library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_fitness", {

    if (length(system("which blastp", intern = TRUE)) != 1){
     testthat::skip("Skipping garnish_fitness because ncbiblast+ is not in PATH")
     }

    dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_output.txt") %>%
            leepR::garnish_fitness

    testthat::expect_true(dt[!is.na(A)] %>% nrow == 32)
    testthat::expect_true(dt[, nmer %>% unique %>% length] == 308)
    testthat::expect_true(dt[, NeoantigenRecognitionPotential %>% na.omit %>%
                            mean %>% signif(digits = 3)] == 0.167)

    if (file.exists("Output/")) unlink("Output", recursive = TRUE, force = TRUE)

    })
