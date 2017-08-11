library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("run_netMHC", {

   if (!check_pred_tools() %>% unlist %>% all) {
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
  }

  # load test data
  # download test files
  dt <- structure(
      list(
      type = c("netMHC", "netMHCpan", "netMHCIIpan"),
      allele = c("HLA-A0201", "HLA-A02:01", "DRB1_0301"),
      nmer_l = c(8L, 9L, 15L),
      filename = c("netMHC_906b84c7-7442-4371-a94b-f0f1b1c80107.csv",
                   "netMHCpan_af9b5685-733a-49c9-8f46-2ea18eab32f5.csv",
                   "netMHCIIpan_5cf162b8-577c-44ae-9c72-f90a89f28a94.csv"
      ),
      command = c("netMHC -p -l 8 -a HLA-A0201 -f netMHC_906b84c7-7442-4371-a94b-f0f1b1c80107.csv",
      "netMHCpan -p -l 9 -a HLA-A02:01 -f netMHCpan_af9b5685-733a-49c9-8f46-2ea18eab32f5.csv",
      "netMHCIIpan -inptype 1 -length 15 -a DRB1_0301 -f netMHCIIpan_5cf162b8-577c-44ae-9c72-f90a89f28a94.csv"
      )),
      .Names = c("type", "allele",
                 "nmer_l", "filename", "command"),
      class = c("data.table", "data.frame"),
      row.names = c(NA, -3L
  ))
  utils::download.file("http://get.rech.io/netMHC_906b84c7-7442-4371-a94b-f0f1b1c80107.csv",
                       "netMHC_906b84c7-7442-4371-a94b-f0f1b1c80107.csv")
  utils::download.file("http://get.rech.io/netMHCpan_af9b5685-733a-49c9-8f46-2ea18eab32f5.csv",
                       "netMHCpan_af9b5685-733a-49c9-8f46-2ea18eab32f5.csv")
  utils::download.file("http://get.rech.io/netMHCIIpan_5cf162b8-577c-44ae-9c72-f90a89f28a94.csv",
                       "netMHCIIpan_5cf162b8-577c-44ae-9c72-f90a89f28a94.csv")
  # run test
  dto <- run_netMHC(dt) %>% data.table::rbindlist(fill = TRUE)

   testthat::expect_gt(dto %>% nrow, 100)
   testthat::expect_gt(dto %>% length, 10)
   testthat::expect_gt(dto$`affinity(nM)_netMHC` %>%
                        stats::na.omit %>%
                        as.numeric %>%
                        sum, 1000000)
   testthat::expect_true(c("AQSGTPPT",
                           "AYESSEDC",
                           "CSPRDRFL",
                           "CSPWDRFL",
                           "ENYWRKAY") %chin% dto$icore_netMHC %>% all)

    })