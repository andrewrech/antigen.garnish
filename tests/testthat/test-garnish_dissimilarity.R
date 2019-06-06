testthat::test_that("garnish_dissimilarity", {

  v <- c("SIINFEKL")

  st <- v %>% garnish_dissimilarity(db = "human", aval = 52)

  testthat::expect_equal(st[, dissimilarity] %>% format(digits = 5),
                        c("1"))

})
