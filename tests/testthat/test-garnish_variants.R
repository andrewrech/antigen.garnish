testthat::test_that("garnish_variants", {
  d <- test_data_dir()

  # load test data
  dt <- "antigen.garnish_example.vcf" %>%
    file.path(d, .) %>%

    # run test
    garnish_variants()

  testthat::expect_equal(dt %>% class() %>% .[1], "data.table")
  testthat::expect_equal(dt$cDNA_change, c("c.4988C>T"))


  dt <- "antigen.garnish_hg19anno_example.vcf" %>%
    file.path(d, .) %>%

    # run test
    garnish_variants()

  testthat::expect_equal(dt %>% class() %>% .[1], "data.table")
  testthat::expect_equal(
    dt$ensembl_transcript_id[c(1, 10, 16)],
    c("ENST00000623167.1", "ENST00000254821.1", "ENST00000424564.1")
  )
})
