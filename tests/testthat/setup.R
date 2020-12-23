library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

options(mc.cores = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

# load antigen.garnish test helper functions

skip_pred_tools <- function() {
  if (!antigen.garnish:::check_pred_tools() %>%
    unlist() %>%
    all()) {
    testthat::skip("Skipping antigen.garnish:::run_netMHC because prediction tools are not in PATH")
  }

  if (suppressWarnings(system("which blastp 2> /dev/null", intern = TRUE)) %>%
    length() == 0) {
    testthat::skip("Skipping garnish_affinity README example from VCF because ncbiblast+ is not in PATH")
  }
}

test_data_dir <- function() {
  dir <- system.file(package = "antigen.garnish") %>% file.path(., "extdata/testdata")
  return(dir)
}
