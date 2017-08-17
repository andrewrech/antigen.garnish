library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("antigen.garnish neoepitope prediction", {

# visually inspect sliding peptides

library(magrittr)

  dt <- data.table::data.table(
       pep_base = "Y___*___THIS_IS_________*___A_NICE_TEST_______*__X",
       mutant_loc = c(5, 25, 47, 50),
       pep_type = "test",
       var_uuid = c("middle",
                    "back_truncate",
                    "front_truncate",
                    "end")) %>%

    get_nmers %T>% print

    testthat::succeed()
    })
