if(requireNamespace("testthat", quietly = TRUE)){
    library(testthat)
    library(antigen.garnish)
    # debug(test_file)
    # debug(testthat:::test_dir)
    test_check(package = "antigen.garnish",
               wrap = FALSE,
               stop_on_failure = FALSE)
}

