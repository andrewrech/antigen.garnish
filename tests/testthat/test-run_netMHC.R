testthat::test_that("run_netMHC", {

  list.files(pattern = "netMHC.*csv") %>% file.remove
  on.exit(list.files(pattern = "netMHC_.*csv") %>% file.remove)

  skip_pred_tools()

  # load test data
    dt <- structure(
        list(
        type = c("netMHC", "netMHCpan", "netMHCIIpan"),
        allele = c("HLA-A0201", "HLA-A02:01", "DRB1_0301"),
        nmer_l = c(8L, 9L, 15L),
        filename = c("antigen.garnish_example_netMHC_906b84c7-7442-4371.csv",
                     "antigen.garnish_example_netMHCIIpan_af9b5685-733a-49c9.csv",
                     "antigen.garnish_example_netMHCpan_5cf162b8-577c-44ae.csv"
        ),
        command = c("netMHC -p -l 8 -a HLA-A0201 -f antigen.garnish_example_netMHC_906b84c7-7442-4371.csv",
        "netMHCpan -p -l 9 -a HLA-A02:01 -f antigen.garnish_example_netMHCIIpan_af9b5685-733a-49c9.csv",
        "netMHCIIpan -inptype 1 -length 15 -a DRB1_0301 -f antigen.garnish_example_netMHCpan_5cf162b8-577c-44ae.csv"
        )),
        .Names = c("type", "allele",
                   "nmer_l", "filename", "command"),
        class = c("data.table", "data.frame"),
        row.names = c(NA, -3L
    ))
    utils::download.file("http://get.rech.io/antigen.garnish_example_netMHC_906b84c7-7442-4371.csv",
                         "antigen.garnish_example_netMHC_906b84c7-7442-4371.csv")
    utils::download.file("http://get.rech.io/antigen.garnish_example_netMHCIIpan_af9b5685-733a-49c9.csv",
                         "antigen.garnish_example_netMHCIIpan_af9b5685-733a-49c9.csv")
    utils::download.file("http://get.rech.io/antigen.garnish_example_netMHCpan_5cf162b8-577c-44ae.csv",
                         "antigen.garnish_example_netMHCpan_5cf162b8-577c-44ae.csv")
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