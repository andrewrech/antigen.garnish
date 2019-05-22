testthat::test_that("iedb_score", {

  v <- c("SIINFEKL", "ILAKFLHWL", "ILRGSVAHK")

  ie <- v %>% iedb_score(db = "mouse")

  testthat::expect_equal(ie[, iedb_score] %>% format(digits = 5),
                        c("1.0000e+00", "9.0884e-07", "1.0000e+00"))

})
