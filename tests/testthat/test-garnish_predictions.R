library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_predictions", {

   if (!check_pred_tools() %>% unlist %>% all) {
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
    }

# copy of README

library(magrittr)

    # download an example VCF
    dt <- "antigen.garnish_example.vcf" %T>%
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
    garnish_summary

    # does antigen.garnish work?
    testthat::expect_equal(dt$variants, 3)
    testthat::expect_equal(dt$alt_neos, 2)
    testthat::expect_gt(dt$predictions, 230)

    if (file.exists("antigen.garnish_example.vcf"))
    file.remove("antigen.garnish_example.vcf")
    })
