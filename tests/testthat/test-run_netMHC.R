testthat::test_that("run_netMHC", {

  list.files(pattern = "ex_netMHC.*csv") %>% file.remove
  on.exit(list.files(pattern = "netMHC_.*csv") %>% file.remove)

  skip_pred_tools()
  d <- test_data_dir()

  # load test data
    dt <- structure(
        list(
        type = c("netMHC", "netMHCpan", "netMHCIIpan"),
        allele = c("HLA-A0201", "HLA-A02:01", "DRB1_0301"),
        nmer_l = c(8L, 9L, 15L),
        filename = c("a.g_ex_netMHC_906b84c7-7442-4371.csv",
                     "a.g_ex_netMHCIIpan_af9b5685-733a-49c9.csv",
                     "a.g_ex_netMHCpan_5cf162b8-577c-44ae.csv"
        ),
        command = c("netMHC -p -l 8 -a HLA-A0201 -f a.g_ex_netMHC_906b84c7-7442-4371.csv",
        "netMHCpan -p -l 9 -a HLA-A02:01 -f a.g_ex_netMHCIIpan_af9b5685-733a-49c9.csv",
        "netMHCIIpan -inptype 1 -length 15 -a DRB1_0301 -f a.g_ex_netMHCpan_5cf162b8-577c-44ae.csv"
        )),
        .Names = c("type", "allele",
                   "nmer_l", "filename", "command"),
        class = c("data.table", "data.frame"),
        row.names = c(NA, -3L
    )) %>%
    # necessary to prevent invalid internal selfref warning with set syntax
    data.table::setDT()

    input <- file.path(d, c("a.g_ex_netMHC_906b84c7-7442-4371.csv",
    "a.g_ex_netMHCIIpan_af9b5685-733a-49c9.csv",
    "a.g_ex_netMHCpan_5cf162b8-577c-44ae.csv"
    ))

    # copy to wd
    file.copy(input, basename(input))

  # run test
    dto <- run_netMHC(dt) %>% data.table::rbindlist(fill = TRUE)

   testthat::expect_gt(dto %>% nrow, 100)
   testthat::expect_gt(dto %>% length, 10)
   testthat::expect_gt(dto$`affinity(nM)_netMHC` %>%
                        stats::na.omit() %>%
                        as.numeric %>%
                        sum, 1000000)
   testthat::expect_true(c("AQSGTPPT",
                           "AYESSEDC",
                           "CSPRDRFL",
                           "CSPWDRFL",
                           "ENYWRKAY") %chin% dto$icore_netMHC %>% all)

    })
