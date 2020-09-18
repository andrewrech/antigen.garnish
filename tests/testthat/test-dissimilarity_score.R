testthat::test_that("dissimilarity_score", {

  v <- c("SIINFEKL")

  st <- v %>% dissimilarity_score(db = "human", aval = 52)

  testthat::expect_equal(st[, dissimilarity] %>% format(digits = 5),
                        c("1"))

})
