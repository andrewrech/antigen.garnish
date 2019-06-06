testthat::test_that("SW_align", {

  v1 <- c("SIINFEKL", "SYFPEITHI")
  v2 <- c("SIIPFEKL", "SYFFPEITHI")

  ot <- SW_align(v1, v2)

  testthat::expect_equal(ot, c(30, 42))

})
