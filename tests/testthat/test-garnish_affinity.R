
README <- function() {
  testthat::test_that("garnish_affinity strict README example", {
    skip_pred_tools()

    # extact copy of README

    # load an example VCF
    dir <- system.file(package = "antigen.garnish") %>%
      file.path(., "extdata/testdata")

    dt <- "antigen.garnish_example.vcf" %>%
      file.path(dir, .) %>%

      # extract variants
      garnish_variants(.) %>%

      # add space separated MHC types
      # see list_mhc() for nomenclature of supported alleles

      .[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")] %>%

      # predict neoantigens
      garnish_affinity(.)

    # this length can vary if prediction tool versions drop allele coverage
    # eg newest netMHC drops A01:47, so don't check row numbers of table or
    # number of unique MHC types
    testthat::expect_true(dt[
      pep_type == "mutnfs",
      nmer %>% unique() %>% length()
    ] == 91)
    testthat::expect_true(dt[, nmer %>% unique() %>% length()] == 184)
  })
}

transcripts <- function() {
  testthat::test_that("garnish_affinity from transcripts, diverse MHC", {
    skip_pred_tools()

    # load test data
    dt <- data.table::data.table(
      sample_id = "test",
      ensembl_transcript_id =
        c("ENSMUST00000128119.1"),
      cDNA_change = c("c.4988C>T"),
      MHC = c("HLA-A*02:01 HLA-DRB1*14:67")
    ) %>%
      # run test
      garnish_affinity(
        blast = FALSE,
        save = FALSE
      )

    testthat::expect_equal(
      dt[, .N, by = "MHC"] %>% .[order(MHC)],
      structure(list(
        MHC = c("HLA-A*02:01", "HLA-DRB1*14:67"),
        N = c(153L, 30L)
      ),
      row.names = c(NA, -2L),
      class = c("data.table", "data.frame")
      )
    )
  })
}

peptides <- function() {
  testthat::test_that("garnish_affinity assemble from peptides", {

    # load test data
    dt <- data.table::data.table(
      sample_id = "test",
      pep_mut = "AAVMILKWTFRAKINGDEHKRDEF",
      mutant_index = c("7 13 14", "all", NA),
      MHC = c("H-2-Kb")
    )
    # run test data
    dto <- garnish_affinity(dt,
      predict = FALSE,
      blast = FALSE
    )

    testthat::expect_equal(
      dto$nmer %>% unique() %>% length(),
      109
    )
  })
}

peptides_wt <- function() {
  testthat::test_that("garnish_affinity assemble from peptides with wild-type", {

    # load test data
    skip_pred_tools()

    d <- test_data_dir()

    # load test data
    dt <- file.path(d, "antigen.garnish_example_peptide_with_WT_input.txt") %>%
      data.table::fread()

    w <- try(dt[sample_id %like% "err.*consecutive"] %>%
      garnish_affinity(blast = FALSE, predict = FALSE))

    testthat::expect_equal(class(w), "try-error")

    w <- try(dt[sample_id == "err_pep_mut"] %>%
      garnish_affinity(blast = FALSE, predict = FALSE))

    testthat::expect_equal(class(w), "try-error")

    w <- try(dt[sample_id == "err_pep_wt"] %>%
      garnish_affinity(blast = FALSE, predict = FALSE))

    testthat::expect_equal(class(w), "try-error")

    # dropping stop gained on dual peptide input throws warning, suppress
    w <- suppressWarnings(
      try(dt[sample_id %like% "stop_gained"] %>%
        garnish_affinity(blast = FALSE, predict = FALSE))
    )

    # expect "no variants for peptide generation" character here
    testthat::expect_equal(class(w), "character")

    # run test data
    dto <- garnish_affinity(dt[!sample_id %like% "^err|^stop"], blast = FALSE, predict = FALSE)

    a <- dto[!is.na(dai_uuid) & pep_type == "wt",
      nmer %>% unique() %>% length(),
      by = "sample_id"
    ]

    b <- dto[!is.na(dai_uuid) & pep_type != "wt",
      nmer %>% unique() %>% length(),
      by = "sample_id"
    ]

    c <- merge(a, b, by = "sample_id")

    testthat::expect_equal(c[, V1.x], c[, V1.y])
    testthat::expect_equal(dto[, sample_id %>% unique() %>% length()], 2)
  })
}

RNA_test <- function() {
  testthat::test_that("garnish_affinity from transcripts, with RNA", {
    skip_pred_tools()

    # load test data
    dt <- data.table::data.table(
      sample_id = "test",
      ensembl_transcript_id =
        c(
          "ENSMUST00000128119.1",
          "ENSMUST00000044250.1",
          "ENSMUST00000018743.1"
        ),
      cDNA_change = c(
        "c.4988C>T",
        "c.1114T>G",
        "c.718T>A"
      ),
      MHC = "H-2-Kb"
    )

    data.table::data.table(
      id =
        c(
          "ENSMUST00000128119.1",
          "ENSMUST00000044250.1",
          "ENSMUST00000018743.1"
        ),
      test = c(100, 5, 1)
    ) %>%
      data.table::fwrite("antigen.garnish_rna_temp.txt",
        sep = "\t",
        quote = FALSE, row.names = FALSE
      )

    # run test
    dt <- garnish_affinity(dt,
      counts = "antigen.garnish_rna_temp.txt",
      min_counts = 10,
      blast = FALSE,
      save = FALSE
    )

    testthat::expect_true(dt$ensembl_transcript_id %>%
      unique() == "ENSMUST00000128119.1")
  })
}

README()
transcripts()
peptides()
peptides_wt()
RNA_test()

list.dirs(path = ".") %>%
  stringr::str_extract("ag_[a-f0-9]{18}") %>%
  na.omit() %>%
  lapply(., function(dir) {
    unlink(dir, recursive = TRUE, force = TRUE)
  })
