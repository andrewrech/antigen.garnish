testthat::test_that("antigen.garnish:::make_sw_alignment", {

  v1 <- c("SIINFEKL", "SYFPEITHI")
  v2 <- c("SIIPFEKL", "SYFFPEITHI")

  ot <- antigen.garnish:::make_sw_alignment(v1, v2)

  testthat::expect_equal(ot, c(30, 42))

})
