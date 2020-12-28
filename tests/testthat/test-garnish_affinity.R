
README <- function() {
  testthat::test_that("garnish_affinity README examples", {
    skip_pred_tools()

    # load an example VCF
    dir <- system.file(package = "antigen.garnish") %>%
      file.path(., "extdata/testdata")

    file <- file.path(dir, "TUMOR.vcf")

    # extract variants
    dt <- garnish_variants(file)

    # add space separated MHC types
    # see antigen.garnish:::list_mhc() for nomenclature of supported alleles
    # MHC may also be set to "all_human" or "all_mouse" to use all supported alleles

    dt[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")]

    # predict neoantigens
    dt %<>% garnish_affinity(.)

    # check overall peptide creation
    testthat::expect_equal(
      dt[, .N, by = pep_type],
      structure(list(pep_type = c("wt", "mutnfs", "mut_other"), N = c(
        1572L,
        1212L,
        2079L
      )), row.names = c(NA, -3L), class = c(
        "data.table",
        "data.frame"
      )) %>% data.table::as.data.table(.)
    )

    # check unique nmers
    testthat::expect_true(dt[, nmer %>% unique() %>% length()] == 401)

    # check that peptides are created correctly
    pep_check <- dt[, .SD,
      .SDcols = c("cDNA_change", "pep_wt", "pep_mut")
    ] %>% unique()

    # check missense
    check_mis_wt <- pep_check[
      cDNA_change == "c.5188C>T",
      strsplit(pep_wt, split = "")
    ]$V1
    check_mis_mut <- pep_check[
      cDNA_change == "c.5188C>T",
      strsplit(pep_mut, split = "")
    ]$V1

    dif <- vector() %>% as.numeric(.)
    for (n in 1:length(check_mis_mut)) {
      if (check_mis_mut[n] != check_mis_wt[n]) {
        dif %<>% c(n)
      }
    }

    testthat::expect_equal(dif, 1730)

    # check deletion
    pep_check <- dt[, .SD, .SDcols = c("cDNA_change", "pep_wt", "pep_mut")] %>% unique()

    check_del_wt <- pep_check[
      cDNA_change == "c.3327_3347delTGCTGCAGCTGCAGCTGCAGC",
      strsplit(pep_wt, split = "")
    ]$V1
    check_del_mut <- pep_check[
      cDNA_change == "c.3327_3347delTGCTGCAGCTGCAGCTGCAGC",
      strsplit(pep_mut, split = "")
    ]$V1

    dif <- vector() %>% as.numeric(.)
    for (n in 1:length(check_del_mut)) {
      if (check_del_mut[n] != check_del_wt[n]) {
        dif %<>% c(n)
      }
    }

    testthat::expect_equal(dif %>% .[1], 1113)

    # check fs
    check_fs_wt <- pep_check[
      cDNA_change == "c.5728_5729insAC",
      strsplit(pep_wt, split = "")
    ]$V1
    check_fs_mut <- pep_check[
      cDNA_change == "c.5728_5729insAC",
      strsplit(pep_mut, split = "")
    ]$V1

    dif <- vector() %>% as.numeric(.)
    for (n in 1:length(check_fs_mut)) {
      if (check_fs_mut[n] != check_fs_wt[n]) {
        dif %<>% c(n)
      }
    }

    testthat::expect_equal(dif, c(
      1910, 1911, 1912, 1913, 1914, 1915, 1916, 1917, 1918, 1919,
      1920, 1921, 1922
    ))
  })

  v <- c("SIINFEKL", "ILAKFLHWL", "GILGFVFTL")

  # calculate foreignness score
  iedb <- v %>% foreignness_score(db = "human")

  testthat::expect_equal(iedb$foreignness_score, 1)

  # calculate dissimilarity
  dis <- v %>% dissimilarity_score(db = "human")

  testthat::expect_equivalent(dis$dissimilarity, c(0.16455097200822866, 0, 2.6688207199754288e-11))
}

transcripts <- function() {
  testthat::test_that("garnish_affinity from transcripts, diverse MHC", {
    skip_pred_tools()

    # load test data
    dt <- data.table::data.table(
      sample_id = "test",
      transcript_id =
        c("ENSMUST00000128119.1"),
      cDNA_change = c("c.4988C>T"),
      MHC = c("HLA-A*02:01 HLA-E*01:03", "HLA-DQA10402-DQB10511")
    ) %>%
      # run test
      garnish_affinity(
        blast = FALSE,
        save = FALSE
      )

    testthat::expect_true(
      (dt$MHC %chin% c("HLA-A*02:01", "HLA-E*01:03", "HLA-DQA10402-DQB10511")) %>% all()
    )
  })
}

peptides <- function() {
  testthat::test_that("garnish_affinity assemble from peptides", {
    skip_pred_tools()

    # load test data
    dt <- data.table::data.table(
      sample_id = "test",
      pep_mut = "AAVMILKWTFRAKINGDEHKRDEF",
      mutant_index = c("7 13 14", "all", NA),
      MHC = c("H-2-Kb")
    )
    # run test data
    dto <- garnish_affinity(dt,
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
    skip_pred_tools()

    # load test data
    d <- test_data_dir()

    # load test data
    dt <- file.path(d, "antigen.garnish_example_peptide_with_WT_input.txt") %>%
      data.table::fread()

    w <- try(dt[sample_id %like% "err.*consecutive"] %>%
      garnish_affinity(blast = FALSE))

    testthat::expect_equal(class(w), "try-error")

    w <- try(dt[sample_id == "err_pep_mut"] %>%
      garnish_affinity(blast = FALSE))

    testthat::expect_equal(class(w), "try-error")

    w <- try(dt[sample_id == "err_pep_wt"] %>%
      garnish_affinity(blast = FALSE))

    testthat::expect_equal(class(w), "try-error")

    # dropping stop gained on dual peptide input throws warning, suppress
    w <- suppressWarnings(
      try(dt[sample_id %like% "stop_gained"] %>%
        garnish_affinity(blast = FALSE))
    )

    # expect "no variants for peptide generation" character here
    testthat::expect_equal(class(w), "character")

    # run test data
    dto <- garnish_affinity(dt[!sample_id %like% "^err|^stop"], blast = FALSE)

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
      transcript_id =
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
    dt <- try(garnish_affinity(dt,
      counts = "antigen.garnish_rna_temp.txt",
      min_counts = 1000,
      blast = FALSE,
      save = FALSE
    ), silent = TRUE)

    testthat::expect_equal(dt %>% class(), "try-error")
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
