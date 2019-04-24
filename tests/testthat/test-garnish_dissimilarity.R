testthat::test_that("garnish_dissimilarity", {

  v <- c("SIINFEKL")

  st <- v %>% garnish_dissimilarity(db = "human")

  testthat::expect_equal(st[, dissimilarity] %>% is.na,
                        TRUE)

})
