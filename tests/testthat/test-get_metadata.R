testthat::test_that("get_metadata", {

  # load test data
  dto <- data.table::data.table(
    transcript_id = c(
      "ENST00000256078.10",
      "ENSMUST00000109987.1"
    )
  ) %>%
    # run tets
    get_metadata()

  testthat::expect_equal(
    dto %>% names(),
    c(
      "transcript_id",
      "peptide",
      "codingTo",
      "codingFrom",
      "cDNA",
      "coding",
      "pepFromCoding",
      "description",
      "chromosome_name",
      "ensembl_gene_id"
    )
  )

  testthat::expect_equal(dto[transcript_id == "ENST00000256078.10", peptide], "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGCVKIKKCIIM*")
})
