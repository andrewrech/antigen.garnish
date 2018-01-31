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
        testthat::expect_true(dt[!is.na(R), nmer %>% unique %>% length] == 252)
        testthat::expect_true(dt %>% nrow == 308)
        testthat::expect_true(dt[, nmer %>% unique %>% length] == 308)
        testthat::expect_true(dt[, (DAI * R) %>% mean(na.rm = TRUE) %>% signif(digits = 3)] == 0.583)
        testthat::expect_true(dt[nchar(nmer) == 9, (DAI * R) %>% mean(na.rm = TRUE) %>% signif(digits = 3)] == 0.610)

    })
