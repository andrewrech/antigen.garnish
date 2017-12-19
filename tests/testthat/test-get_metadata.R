testthat::test_that("get_metadata", {

  # load test data
    dto <- data.table::data.table(
             ensembl_transcript_id = c("ENSMUST00000128119",
                                        "ENSMUST00000044250",
                                        "ENST00000256078")) %>%
  # run tets
   get_metadata

 testthat::expect_equal(dto %>% names,
        c("ensembl_transcript_id",
          "external_gene_name",
          "ensembl_gene_id",
          "description",
          "chromosome_name",
          "start_position",
          "end_position",
          "transcript_start",
          "transcript_end",
          "coding",
          "peptide")
)

  # load test data
    dto <- data.table::data.table(
             ensembl_transcript_id = c("ENSMUST00000128119",
                                        "ENSMUST00000044250",
                                        "ENST00000256078")) %>%
  # run tets
   get_metadata(humandb = "GRCh37",
                mousedb = "GRCm37")

 testthat::expect_equal(dto %>% names,
        c("ensembl_transcript_id",
          "external_gene_name",
          "ensembl_gene_id",
          "description",
          "chromosome_name",
          "start_position",
          "end_position",
          "transcript_start",
          "transcript_end",
          "coding",
          "peptide")
)

    })
