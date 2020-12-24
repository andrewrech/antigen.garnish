library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

testthat::test_that("antigen.garnish neoantigen prediction", {

# visually inspect sliding peptides

  # generate a fake peptide
    dt <- data.table::data.table(
       pep_base = "Y___*___THIS_IS_________*___A_CODE_TEST!______*__X",
       mutant_index = c(5, 25, 47, 50),
       pep_type = "test",
       var_uuid = c(
                    "front_truncate",
                    "middle",
                    "back_truncate",
                    "end")) %>%
  # create nmers
    make_nmers %T>% print

    testthat::succeed()
    })
