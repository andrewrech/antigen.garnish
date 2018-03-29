if(requireNamespace("testthat", quietly = TRUE)){
    library(testthat)
    library(antigen.garnish)
    test_check(package = "antigen.garnish",
               wrap = FALSE,
               stop_on_failure = FALSE)
}

