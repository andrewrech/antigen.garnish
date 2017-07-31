library(testthat)
library(antigen.garnish)
library(data.table)

test_check("antigen.garnish")

test_that("antigen.garnish pipeline", {

mhc_dt <- data.table::data.table(
            transcript_affected = c("ENST00000256078", "ENST00000256078"),
            sample_id = c("test_sample_1", "test_sample_2"),
            aa_mutation = c("G12D", "G13D"),
            MHC = c("HLA-A*02:01 HLA-A*03:01 HLA-DRB10301 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67"))

mhc_dt <- garnish_predictions(mhc_dt)

dt <- garnish_summary(mhc_dt)

testthat::compare(dt,
  structure(list(
                 sample_id = c("test_sample_1", "test_sample_2"),
                 priority_neos = c(0L, 0L),
                 classic_neos = c(0L, 0L),
                 alt_neos = c(0L, 0L),
                 alt_neos_top = c(5.20144265891428, 1.56599284273263),
                 classic_neos_top = c(0.00954700205453025, 0.00954700205453025),
                 binders = c(11L, 9L),
                 peptides = c(74L, 76L),
                 predictions = c(148L, 152L)),
                   .Names = c("sample_id", "classic_neos", "alt_neos", "alt_neos_top",
                   "classic_neos_top", "binders", "peptides", "predictions"),
                   row.names = c(NA, -2L), class = c("data.table", "data.frame"))
          )})
