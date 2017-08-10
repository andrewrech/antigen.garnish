library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("antigen.garnish neoepitope prediction", {

# copy of README

library(magrittr)

    # download an example VCF
    "antigen.garnish_example.vcf" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

    # extract variants
    garnish_variants %>%

    # add MHC types
        .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-DRB1*14:67",
                    "HLA-A*03:01 HLA-DRB1*03:01")] %>%

    # predict neoepitopes
    garnish_predictions %>%

    # summarize predictions
    garnish_summary %T>%

    print %>%

    # does antigen.garnish work?
    testthat::compare(
       data.table::fread("http://get.rech.io/antigen.garnish_example_summary.txt"))
    })
