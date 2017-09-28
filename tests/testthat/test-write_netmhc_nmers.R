library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("write_netmhc_nmers", {

  list.files(pattern = "netMHC.*csv") %>% file.remove
  on.exit(list.files(pattern = "netMHC_.*csv") %>% file.remove)

  # load test data
  dt <- data.table::data.table(netMHC = c("HLA-A0201", "HLA-A0201", "HLA-A0201",
  "HLA-A0201", "HLA-A0201", "HLA-A0201", "HLA-A0201", "HLA-A0201",
  "HLA-A0201", "HLA-A0201"), nmer = c("AQSGTPPT", "AYESSEDC", "DDENYWRK",
  "ENYWRDDD", "GAQSGTPP", "GTWVSGAQ", "GTWVSGVQ", "GVQSGTPP", "KAYESSED",
  "KSYESSED"), nmer_l = c(8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L
  ))

  # run test
  dto <-  write_netmhc_nmers(dt, "netMHC")

  testthat::expect_true(
    data.table::fread(dto$filename,
          header = FALSE)$V1 %chin%
    c("AQSGTPPT", "AYESSEDC", "DDENYWRK",
    "ENYWRDDD", "GAQSGTPP", "GTWVSGAQ",
    "GTWVSGVQ", "GVQSGTPP", "KAYESSED",
    "KSYESSED") %>% all
    )

  testthat::expect_equal(dto$allele %>% unique, "HLA-A0201")
    })
