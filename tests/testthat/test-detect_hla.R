library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

testthat::test_that("write_nmers", {

  # load test data
  alleles <- data.table::rbindlist(
                         list(
  system.file("extdata",
      "netMHC_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHC"],
  system.file("extdata",
      "netMHCpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCpan"],
  system.file("extdata",
      "mhcflurry_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcflurry"],
  system.file("extdata",
      "netMHCII_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCII"],
  system.file("extdata",
      "netMHCIIpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCIIpan"]
                      ))
  dt <- data.table::data.table(netMHCIIpan =
         c("A0201", "A0301", "DRB1_0301", "DRB1_1467"))

  # run test
  dt[, netMHCIIpan := detect_hla(netMHCIIpan, alleles)]$netMHCIIpan %>%

   testthat::expect_equal(c(NA, NA, "DRB1_0301", "DRB1_1467"))
    })
