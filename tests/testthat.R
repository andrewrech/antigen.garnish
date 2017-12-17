if(requireNamespace("testthat", quietly = TRUE)){
    library(testthat)
    library(antigen.garnish)
    test_check("antigen.garnish", wrap = FALSE, stop_on_failure = FALSE)
}