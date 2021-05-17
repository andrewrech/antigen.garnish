testthat::test_that("get_pred_commands", {
  list.files(pattern = "(netMHC|flurry).*-.*-.*\\.csv") %>% file.remove()
  on.exit(list.files(pattern = "(netMHC|flurry).*-.*-.*\\.csv") %>% file.remove())

  # load test data

  dt <- data.table::data.table(
    nmer = c(
      "AQSGTPPT", "AQSGTPPTG", "AQSGTPPTGL",
      "AQSGTPPTGLY", "AQSGTPPTGLYG", "AQSGTPPTGLYGH"
    ),
    MHC = c(
      "HLA-A*02:01 HLA-DRB1*14:67",
      "HLA-A*02:01 HLA-DRB1*14:67",
      "HLA-A*02:01 HLA-DRB1*14:67",
      "HLA-A*02:01 HLA-DRB1*14:67",
      "HLA-A*02:01 HLA-DRB1*14:67",
      "HLA-A*02:01 HLA-DRB1*14:67"
    ),
    nmer_l = 8:13
  )

  # run test

  dto <- dt %>%
    antigen.garnish:::get_pred_commands() %>%
    .[[2]]

  testthat::expect_equal(dto %>% nrow(), 12)
  testthat::expect_equal(
    dto %>% names(),
    c("type", "allele", "nmer_l", "filename", "command")
  )
})
