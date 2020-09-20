testthat::test_that("dissimilarity_score", {

  v <- c("SIINFEKL")

  check_pred_tools()

  st <- v %>% dissimilarity_score(db = "human", aval = 52)

  testthat::expect_equal(st[, dissimilarity] %>% format(digits = 5),
                        c("1"))

})
