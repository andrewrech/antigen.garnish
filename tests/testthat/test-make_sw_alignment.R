testthat::test_that("make_sw_alignment", {

  v1 <- c("SIINFEKL", "SYFPEITHI")
  v2 <- c("SIIPFEKL", "SYFFPEITHI")

  ot <- make_sw_alignment(v1, v2)

  testthat::expect_equal(ot, c(30, 42))

})
