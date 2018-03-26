testthat::test_that("list_mhc", {

  testthat::expect_equal(
    list_mhc() %>% nrow,
    3874)

    })
