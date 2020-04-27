testthat::test_that("garnish_antigens", {
  skip_pred_tools()
  d <- test_data_dir()

  # load test data
  dt <- file.path(d, "antigen.garnish_PureCN_example_output.txt") %>%
    data.table::fread() %>%
    garnish_antigens()

  # run test

  testthat::expect_equal(dt, structure(
    list(
      sample_id = c("antigen.garnish_test_pureCN.vcf", "antigen.garnish_test_pureCN.vcf", "antigen.garnish_test_pureCN.vcf"),
      nmer = c("KLVVVGAVGV", "LICPYMEPI", "QQLEAAAGV"),
      MHC = c("HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01"),
      protein_change = c("p.Gly12Val", "p.Arg769Cys", "p.Pro1426Leu"),
      cDNA_change = c("c.35G>T", "c.2305C>T", "c.4277C>T"),
      external_gene_name = c("KRAS", "MTOR", "MTOR"),
      Ensemble_score = c(323.874201423834, 132.615664199976, 340.278748742299),
      dissimilarity = c(NA_real_, NA_real_, NA_real_),
      iedb_score = c(NA, 0.999999999973312, NA),
      min_DAI = c(3.23053559356532, 12.505756955261, 3.0504938492974),
      clone_id = c(1L, 2L, 2L),
      cl_proportion = c(0.583196682780358, 0.277452002724779, 0.277452002724779)
    ),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -3L)
  ))
})
