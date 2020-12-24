testthat::test_that("prediction output", {
  convert.aa <- function(x) {
    sn <- names(Biostrings::AMINO_ACID_CODE)
    x <- Biostrings::AMINO_ACID_CODE[which(sn == x)]
    unlist(x)
  }

  d <- test_data_dir()


  # load test data
  dt <- data.table::fread(file.path(d, "antigen.garnish_example_affinity_output.txt"))

  # test if nmers are contained in peptides

  dt[pep_type != "wt", test := stringr::str_extract(pattern = nmer, string = pep_mut)]
  dt[pep_type == "wt", test := stringr::str_extract(pattern = nmer, string = pep_wt)]
  dt[, validate := test == nmer,
    by = 1:nrow(dt)
  ]

  testthat::expect_true(dt$validate %>% all())

  # test that DAI peptides are only 1 AA apart

  DAIdt <- dt[!is.na(dt$DAI)][pep_type != "wt"]

  DAIdt[, n.aa.mismatch := lapply(1:nrow(DAIdt), function(i) {
    x <- pep_mut[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    y <- pep_wt[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()

    return(sum(x != y))
  }) %>% unlist()]

  testthat::expect_true(all(DAIdt$n.aa.mismatch == 1))

  # match mutant_index to SnpEff protein change call for all mutants

  dt[, prot.change.from.se := stringr::str_extract_all(string = protein_change, pattern = "[0-9]+") %>% unlist()]

  testthat::expect_true(dt[, mutant_index == prot.change.from.se] %>% all())

  # test that pep_mut matches SnpEff protein change call for missense

  DAIdt[, new_aa_pep := lapply(1:nrow(DAIdt), function(i) {
    x <- pep_mut[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    y <- pep_wt[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    nAA <- which(x != y)

    return(x[nAA])
  }) %>% unlist()]

  DAIdt[, old_aa_pep := lapply(1:nrow(DAIdt), function(i) {
    x <- pep_mut[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    y <- pep_wt[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    nAA <- which(x != y)

    return(y[nAA])
  }) %>% unlist()]

  DAIdt[, old_aa_pep := convert.aa(old_aa_pep),
    by = 1:nrow(DAIdt)
  ][, new_aa_pep := convert.aa(new_aa_pep),
    by = 1:nrow(DAIdt)
  ]

  DAIdt[, old_aa_se :=
    stringr::str_extract(
      string = protein_change,
      pattern = "(?<=(p\\.))[A-Z][a-z]{2}(?=[0-9])"
    ) %>% unlist()]
  DAIdt[, new_aa_se :=
    stringr::str_extract(
      string = protein_change,
      pattern = "(?<=[0-9])[A-Z][a-z]{2}$"
    ) %>% unlist()]

  testthat::expect_true(DAIdt[, old_aa_pep == old_aa_se] %>% all())
  testthat::expect_true(DAIdt[, new_aa_pep == new_aa_se] %>% all())

  # check wt and mut peptides are registered correctly

  # getting missense AA change

  DAIdt[, new_aa_pep := lapply(1:nrow(DAIdt), function(i) {
    x <- pep_mut[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    y <- pep_wt[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    nAA <- which(x != y)

    return(x[nAA])
  }) %>% unlist()]

  DAIdt[, old_aa_pep := lapply(1:nrow(DAIdt), function(i) {
    x <- pep_mut[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    y <- pep_wt[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    nAA <- which(x != y)

    return(y[nAA])
  }) %>% unlist()]

  # deriving wt from mutant peptide

  DAIdt[, pep_wt2 := lapply(1:nrow(DAIdt), function(i) {
    x <- pep_mut[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    if (x[mutant_index[i]] == new_aa_pep[i]) {
      x[mutant_index[i]] <- old_aa_pep[i]
    } else {
      return(NA)
    }

    return(paste(x, collapse = ""))
  }) %>% unlist()]

  testthat::expect_true(DAIdt[, pep_wt == pep_wt2] %>% all())

  # deriving mutant from wt peptide

  DAIdt[, pep_mut2 := lapply(1:nrow(DAIdt), function(i) {
    x <- pep_wt[i] %>%
      strsplit(
        split = "",
        fixed = TRUE
      ) %>%
      unlist()
    if (x[mutant_index[i]] == old_aa_pep[i]) {
      x[mutant_index[i]] <- new_aa_pep[i]
    } else {
      return(NA)
    }

    return(paste(x, collapse = ""))
  }) %>% unlist()]

  testthat::expect_true(DAIdt[, pep_mut == pep_mut2] %>% all())

  # compare actual to expected number of peptides generated

  ndt <- DAIdt[, c("nmer_uuid", "var_uuid")] %>% unique()

  if (any(ndt[, .N,
    by = var_uuid
  ]$N < 77)) {
    ndt[, less := .N < 77,
      by = var_uuid
    ]
    vect <- ndt[less == TRUE]$var_uuid
    test_dt <- DAIdt[var_uuid %chin% vect]
    test_dt[, nchar := nchar(pep_mut)]
    test_dt[, min := min(c(mutant_index, (nchar - mutant_index))),
      by = 1:nrow(test_dt)
    ]

    testthat::expect_true(all(test_dt$min < 14))
  }
})
