testthat::test_that("garnish_antigens", {
  skip_pred_tools()
  d <- test_data_dir()

  # load test data
  dt <- file.path(d, "antigen.garnish_example_output2.txt") %>%
    data.table::fread() %>%
    garnish_antigens()

  # run test

  testthat::expect_equal(dt, structure(list(sample_id = c(
    "antigen.garnish_test_pureCN.vcf",
    "antigen.garnish_test_pureCN.vcf", "antigen.garnish_test_pureCN.vcf",
    "antigen.garnish_test_pureCN.vcf", "antigen.garnish_test_pureCN.vcf"
  ), nmer = c(
    "INNKLQQLE", "KLQQLEAAA", "LEAAAGVLE", "LQQLEAAAG",
    "NKLQQLEAA"
  ), MHC = c(
    "HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01",
    "HLA-A*02:01", "HLA-A*02:01"
  ), ensembl_transcript_id = c(
    "ENST00000361445",
    "ENST00000361445", "ENST00000361445", "ENST00000361445", "ENST00000361445"
  ), protein_change = c(
    "p.Pro1426Leu", "p.Pro1426Leu", "p.Pro1426Leu",
    "p.Pro1426Leu", "p.Pro1426Leu"
  ), cDNA_change = c(
    "c.4277C>T",
    "c.4277C>T", "c.4277C>T", "c.4277C>T", "c.4277C>T"
  ), clone_id = c(
    2L,
    2L, 2L, 2L, 2L
  ), Ensemble_score = c(
    33.537419999999997, 12.080912700000001,
    5.0948799999999999, 30.55461, 10.777609999999999
  ), dissimilarity = c(
    0.5,
    NA, 0, 0.999, 0.40000000000000002
  ), iedb_score = c(
    0.66666666699999999,
    NA, 4.5299999999999999e-07, 0, NA
  ), min_DAI = c(
    0.102436914,
    2.5667469110000001, 14.103280010000001, 0.13068237199999999,
    11.13865477
  ), Recognition_Features = c(
    "binding affinity; foreignness; dissimilarity",
    "binding affinity", "binding affinity; foreignness; agretopicity",
    "binding affinity; dissimilarity", "binding affinity; agretopicity; dissimilarity"
  )), row.names = c(NA, -5L), class = c("data.table", "data.frame")))
})
