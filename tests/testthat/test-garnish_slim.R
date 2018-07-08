test_data_dir <- function(){

  dir <- system.file(package = "antigen.garnish") %>% file.path(., "extdata/testdata")

  return(dir)

}


testthat::test_that("garnish_slim", {

  d <- test_data_dir()

  # load test data
  # need to make test file
    dt <- data.table::fread(file.path(d, "antigen.garnish_PureCN_example_output.txt"))

    dto <- garnish_slim(dt)


    class_I_prediction_cols <- c(
        "mhcflurry_prediction", "mhcflurry_prediction_low",
        "mhcflurry_prediction_high", "mhcflurry_prediction_percentile",
        "mhcnuggets_pred_lstm", "mhcnuggets_pred_gru",
        "affinity(nM)_netMHC", "affinity(nM)_netMHCpan",
        "%Rank_netMHC", "%Rank_netMHCpan"
        )

    class_II_prediction_cols <- c(
        "affinity(nM)_netMHCII", "affinity(nM)_netMHCIIpan",
        "%Random_netMHCII",  "%Rank_netMHCIIpan"
        )


  # run test
     testthat::expect_equal(dto$nmers, 154)
     testthat::expect_equal(dto$mhc_binders_class_I, 3)
     testthat::expect_equal(dto$variants, 2)
     testthat::expect_equal(ncol(dto), 30)


     # Run test in peptide input
     dt_pep <- data.table::data.table(
           sample_id = "test",
           pep_mut = "MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP",
           mutant_index = "12",
           MHC = "HLA-A*02:01") %>%
     antigen.garnish::garnish_affinity

     testthat::expect_equal
    })






















