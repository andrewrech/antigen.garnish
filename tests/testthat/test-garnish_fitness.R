testthat::test_that("garnish_fitness", {

  if (suppressWarnings(system('which blastp 2> /dev/null', intern = TRUE)) %>%
          length == 0)
      testthat::skip("Skipping garnish_fitness because ncbiblast+ is not in PATH")

  if (!check_pred_tools() %>% unlist %>% all){
        testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
        }

    # load test data
        dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_output.txt") %>%
                garnish_fitness

    # run test
        testthat::expect_true(dt[!is.na(A)] %>% nrow == 32)
        testthat::expect_true(dt[, nmer %>% unique %>% length] == 308)
        testthat::expect_true(dt[, NeoantigenRecognitionPotential %>% na.omit %>%
                                mean %>% signif(digits = 3)] == 0.167)

    })
