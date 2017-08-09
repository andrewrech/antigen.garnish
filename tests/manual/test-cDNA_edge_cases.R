library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

testthat::test_that("get_cDNA", {

  # load test data

  dt <- data.table::data.table(
      cDNA_type = c(">", ">", "ins", "del", ">", "del"),
      coding = c(LETTERS %>% paste(collapse = "")),
      cDNA_locs = c(7L, 1L, 14L, 4L, 100L, 1L),
      cDNA_locl = c(7L, 1L, 14L, 7L, 100L, 23L),
      cDNA_seq = c("*", "*", "$$$$$", "", "*", ""))

  # run test
  dto <- get_cDNA(dt)

  # visually inspect
  dto[, `:=`(coding_l = coding %>% nchar,
          coding_mut_l = coding_mut %>% nchar,
          ldelta =  (coding_mut %>% nchar) - (coding %>% nchar)
          )] %>% print

  testthat::succeed()
    })