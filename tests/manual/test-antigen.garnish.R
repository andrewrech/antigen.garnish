library(testthat)
library(antigen.garnish)
library(data.table)

test_that("antigen.garnish neoepitope prediction", {

dt <-
    # load an example VCF
    system.file("extdata",
          "antigen.garnish_example.vcf",
          package = "antigen.garnish") %>%

    # extract variants
    garnish_variants(.)


    # add MHC types
    package_test <- dt$antigen.garnish_input %>%
        .[, MHC := c("HLA-A*02:01 HLA-A*03:01 HLA-DRB10301 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67")] %>%

    # predict
    garnish_predictions(.) %>%

    # sumarize
    garnish_summary(.)


    # does antigen.garnish work?

    testthat::compare(package_test,
    structure(list(sample_id = "tumor",
                    priority_neos = 0L,
                    classic_neos = 0L,
                    alt_neos = 2L,
                    alt_neos_top = 32.5590961308976,
                    classic_neos_top = 0.0086649901329808,
                    binders = 7L,
                    peptides = 231L,
                    predictions = 462L),
                    .Names = c("sample_id",
                    "priority_neos", "classic_neos", "alt_neos", "alt_neos_top",
                    "classic_neos_top", "binders", "peptides", "predictions"),
                    row.names = c(NA, -1L),
                    class = c("data.table", "data.frame"))
    )
    })
