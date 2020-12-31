testthat::test_that("collate_netMHC", {

  # load test data
  file.copy(
    file.path(test_data_dir(), "antigen.garnish_collate_example.txt"),
    "antigen.garnish_collate_example.txt"
  )

  esl <- list(
    list(
      "netMHC -p -l 8 -a HLA-A0201 -f netMHC_906b84c7-7442-4371-a94b-f0f1b1c80107.csv",
      "antigen.garnish_collate_example.txt"
    )
  )

  # run test
  dt <- esl %>% antigen.garnish:::collate_netMHC()

  testthat::expect_equal(dt[[1]] %>% nrow(), 48)
  testthat::expect_equal(dt[[1]] %>% length(), 15)
  testthat::expect_true(c(
    "AQSGTPPT",
    "AYESSEDC",
    "CSPRDRFL",
    "CSPWDRFL",
    "ENYWRKAY"
  ) %chin% dt[[1]]$icore_netMHC %>% all())
})
