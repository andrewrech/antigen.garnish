testthat::test_that("garnish_clonality", {

  # load test data
    dt <- "http://get.rech.io/antigen.garnish_pureCN_example_output.txt" %>%
              data.table::fread %>%
              .[, c("clone_prop", "clone_id", "garnish_score") := NULL]
  # run test

  dt %<>% garnish_clonality

  testthat::expect_equal(dto[!is.na(garnish_score), garnish_score %>% unique %>% length], 1)
  testthat::expect_equal(dto[, clone_id %>% as.numeric %>% range(na.rm = TRUE)], c(1,2))
  testthat::expect_equal(dto[!is.na(clone_prop), clone_prop %>% range(na.rm = TRUE)], c(0.2720, 1.0484))

# now check AF instead of CF


dt %>% data.table::setnames("CELLFRACTION", "AF")

dt %<>% clone_wars


testthat::expect_equal(dt[!is.na(garnish_score), garnish_score %>% unique %>% length], 1)
testthat::expect_equal(dt[, clone_id %>% as.numeric %>% range(na.rm = TRUE)], c(1,2))
testthat::expect_equal(dt[!is.na(clone_prop), clone_prop %>% range(na.rm = TRUE)] %>% signif(digits = 3), c(0.259, 1.000))

    })
