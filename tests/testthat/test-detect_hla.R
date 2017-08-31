library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("detect_mhc", {

  # load test data
  alleles <- data.table::rbindlist(
                         list(
  system.file("extdata",
      "netMHC_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHC"],
  system.file("extdata",
      "netMHCpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCpan"],
  system.file("extdata",
      "mhcflurry_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcflurry"],
  system.file("extdata",
      "netMHCII_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCII"],
  system.file("extdata",
      "netMHCIIpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCIIpan"],
  system.file("extdata",
                "mhcnuggets_gru_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcnuggets_gru"],
  system.file("extdata",
              "mhcnuggets_lstm_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcnuggets_lstm"]
))

  dt <- data.table::data.table(netMHCIIpan =
         c("A0201", "A0301", "DRB1_0301", "DRB1_1467"))

  # run test
  dt[, netMHCIIpan := detect_mhc(netMHCIIpan, alleles)]$netMHCIIpan %>%

   testthat::expect_equal(c(NA, NA, "DRB1_0301", "DRB1_1467"))
    })
