library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)

test_that("write_nmers", {

  on.extit({
    list.files(pattern = "netMHC_.*csv") %>%
    file.remove})

  # load test data

  dt <- structure(list(netMHC = c("HLA-A0201", "HLA-A0201", "HLA-A0201",
  "HLA-A0201", "HLA-A0201", "HLA-A0201", "HLA-A0201", "HLA-A0201",
  "HLA-A0201", "HLA-A0201"), nmer = c("AQSGTPPT", "AYESSEDC", "ENYWRKAY",
  "ENYWRKSY", "GAQSGTPP", "GTWVSGAQ", "GTWVSGVQ", "GVQSGTPP", "KAYESSED",
  "KSYESSED"), nmer_l = c(8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L
  )), .Names = c("netMHC", "nmer", "nmer_l"), sorted = c("netMHC",
  "nmer_l"), class = c("data.table", "data.frame"), row.names = c(NA,
  -10L))

  # run test
  dto <-  write_nmers(dt, "netMHC")

  testthat::expect_equal(
    data.table::fread(dto$filename,
          header = FALSE)$V1,
    c("AQSGTPPT", "AYESSEDC", "ENYWRKAY",
      "ENYWRKSY", "GAQSGTPP", "GTWVSGAQ",
      "GTWVSGVQ", "GVQSGTPP", "KAYESSED",
      "KSYESSED"))
  testthat::expect_equal(dto$allele %>% unique, "HLA-A0201")

    })
