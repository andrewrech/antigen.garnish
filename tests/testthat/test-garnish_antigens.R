testthat::test_that("garnish_antigens", {
  skip_pred_tools()
  d <- test_data_dir()

  # load test data
  dt <- file.path(d, "antigen.garnish_PureCN_example_output2.txt") %>%
    data.table::fread() %>%
    garnish_antigens()

  # run test

  testthat::expect_equal(dt, structure(
    list(
      sample_id = c("antigen.garnish_test_pureCN.vcf",
      "antigen.garnish_test_pureCN.vcf", "antigen.garnish_test_pureCN.vcf",
      "antigen.garnish_test_pureCN.vcf", "antigen.garnish_test_pureCN.vcf"),
      nmer = c("INNKLQQLE", "KLQQLEAAA", "LEAAAGVLE", "LQQLEAAAG", "NKLQQLEAA"),
      MHC = c("HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01",
      "HLA-A*02:01"),
      protein_change = c("p.Pro1426Leu", "p.Pro1426Leu", "p.Pro1426Leu",
      "p.Pro1426Leu", "p.Pro1426Leu"),
      cDNA_change = c("c.4277C>T", "c.4277C>T", "c.4277C>T", "c.4277C>T",
      "c.4277C>T"),
      external_gene_name = c("MTOR", "MTOR", "MTOR", "MTOR", "MTOR"),
      clone_id = c(2L, 2L, 2L, 2L, 2L),
      Ensemble_score = c(33.53742, 12.0809127, 5.09488, 30.55461, 10.77761),
      dissimilarity = c(0.5, NA, 0, 0.999, 0.4),
      iedb_score = c(0.666666667, NA, 4.53e-07, 0, NA),
      min_DAI = c(0.102436914, 2.566746911, 14.10328001, 0.130682372,
      11.13865477),
      Recognition_Features = c("binding affinity; foreignness; dissimilarity",
      "binding affinity", "binding affinity; foreignness; agretopicity",
      "binding affinity; dissimilarity",
      "binding affinity; agretopicity; dissimilarity")),
      row.names = c(NA, -5L),
      class = c("data.table", "data.frame")
    )
  )

})
