testthat::test_that("antigen.garnish:::make_cDNA edge cases", {

  # load test data

  dt <- data.table::data.table(
      cDNA_type = c(">", ">", "ins", "del", ">", "del", "delins"),
      coding = c(LETTERS %>% paste(collapse = "")),
      cDNA_locs = c(7L, 1L, 14L, 4L, 100L, 1L, 14L),
      cDNA_locl = c(7L, 1L, 14L, 7L, 100L, 23L, 16L),
      cDNA_seq = c("*", "*", "_____", "", "*", "", "**"))

  # run test
  dto <- antigen.garnish:::make_cDNA(dt)

  # visually inspect
  dto[, .(coding_mut, coding_l = coding %>% nchar,
          coding_mut_l = coding_mut %>% nchar,
          ldelta =  (coding_mut %>% nchar) - (coding %>% nchar)
          )] %>% print

  testthat::expect_equivalent(dto,
    data.table(structure(
      list(
      cDNA_type = c(">", ">", "ins", "del", "del", "delins"),
      coding = c("ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
      cDNA_locs = c(7L,
                    1L,
                    14L,
                    4L,
                    1L,
                    14L),
      cDNA_locl = c(7L,
                    1L,
                    14L,
                    7L,
                    23L,
                    16L),
      cDNA_seq = c("*",
                   "*",
                   "_____",
                   "",
                   "",
                   "**"),
      coding_mut = c("ABCDEF*HIJKLMNOPQRSTUVWXYZ",
                    "*BCDEFGHIJKLMNOPQRSTUVWXYZ",
                     "ABCDEFGHIJKLMN_____OPQRSTUVWXYZ",
                     "ABCHIJKLMNOPQRSTUVWXYZ",
                     "XYZ",
                     "ABCDEFGHIJKLM**QRSTUVWXYZ")),
      .Names = c("cDNA_type",
                 "coding",
                 "cDNA_locs",
                 "cDNA_locl",
                 "cDNA_seq",
                 "coding_mut"),
                 class = c("data.table",
                           "data.frame"))
))
    })
