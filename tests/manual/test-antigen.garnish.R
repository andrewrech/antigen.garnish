library(testthat)
library(antigen.garnish)

test_that("antigen.garnish neoepitope prediction", {

    # download an example VCF
    "antigen.garnish_example.vcf" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

    # extract variants
    antigen.garnish::garnish_variants %>%

    # add MHC types
    .$antigen.garnish_variants %>%
        .[, MHC := c("HLA-A*02:01 HLA-A*03:01 HLA-DRB10301 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67")] %>%

    # predict
    antigen.garnish::garnish_predictions %>%

    # summarize
    antigen.garnish::garnish_summary %T>%

    print %>%

    # does antigen.garnish work?
    testthat::compare(
            data.table::fread("http://get.rech.io/ag_test.csv"))
    })
