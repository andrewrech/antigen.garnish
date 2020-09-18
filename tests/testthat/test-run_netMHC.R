testthat::test_that("run_netMHC", {
  list.files(pattern = "ex_netMHC.*csv") %>% file.remove()
  on.exit(list.files(pattern = "netMHC_.*csv") %>% file.remove())

  skip_pred_tools()
  d <- test_data_dir()

  # load test data
  dt <- structure(
    list(
      type = c("netMHC", "netMHCpan", "netMHCIIpan"),
      allele = c("HLA-A0201", "HLA-A02:01", "DRB1_0301"),
      nmer_l = c(8L, 9L, 15L),
      filename = c(
        "a.g_ex_netMHC_906b84c7-7442-4371.csv",
        "a.g_ex_netMHCIIpan_af9b5685-733a-49c9.csv",
        "a.g_ex_netMHCpan_5cf162b8-577c-44ae.csv"
      ),
      command = c(
        "netMHC -p -l 8 -a HLA-A0201 -f a.g_ex_netMHC_906b84c7-7442-4371.csv",
        "netMHCpan -p -l 9 -a HLA-A02:01 -f a.g_ex_netMHCIIpan_af9b5685-733a-49c9.csv",
        "netMHCIIpan -inptype 1 -length 15 -a DRB1_0301 -f a.g_ex_netMHCpan_5cf162b8-577c-44ae.csv"
      )
    ),
    .Names = c(
      "type", "allele",
      "nmer_l", "filename", "command"
    ),
    class = c("data.table", "data.frame"),
    row.names = c(NA, -3L)
  ) %>%
    # necessary to prevent invalid internal selfref warning with set syntax
    data.table::setDT()

  input <- file.path(d, c(
    "a.g_ex_netMHC_906b84c7-7442-4371.csv",
    "a.g_ex_netMHCIIpan_af9b5685-733a-49c9.csv",
    "a.g_ex_netMHCpan_5cf162b8-577c-44ae.csv"
  ))

  # copy to wd
  file.copy(input, basename(input))

  # run test
  dto <- run_netMHC(dt) %>% data.table::rbindlist(fill = TRUE)

  # this test should fail if any column names change
  # to alert to netMHC dependency changes across versions

  dto_names_test <- c(
    "pos_netMHC", "netMHC", "nmer", "core_netMHC", "Offset_netMHC",
    "I_pos_netMHC", "I_len_netMHC", "D_pos_netMHC", "D_len_netMHC",
    "icore_netMHC", "Identity_netMHC", "1-log50k(aff)_netMHC", "affinity(nM)_netMHC",
    "%Rank_netMHC", "command_netMHC", "pos_netMHCpan", "netMHCpan",
    "core_netMHCpan", "Of_netMHCpan", "Gp_netMHCpan", "Gl_netMHCpan",
    "Ip_netMHCpan", "Il_netMHCpan", "icore_netMHCpan", "Identity_netMHCpan",
    "Score_netMHCpan", "%Rank_netMHCpan", "command_netMHCpan", "pos_netMHCIIpan",
    "netMHCIIpan", "Of_netMHCIIpan", "core_netMHCIIpan", "Core_Rel_netMHCIIpan",
    "Identity_netMHCIIpan", "Score_EL_netMHCIIpan", "%Rank_EL_netMHCIIpan",
    "Exp_Bind_netMHCIIpan", "command_netMHCIIpan"
  )

  testthat::expect_equal(dto$nmer %>% unique() %>% length(), 192)

  dto_mhc <- c("HLA-A02:01", "HLA-A0201", "DRB1_0301")

  dto_mhc_test <- c(
    dto$netMHCpan %>% unique(),
    dto$netMHC %>% unique(),
    dto$netMHCIIpan %>% unique()
  ) %>%
    stats::na.omit(.) %>%
    as.vector(.)

  testthat::expect_equal(dto_mhc, dto_mhc_test)

  dto_Gp <- c("1", "2", "3", "4", "5", "6", "7", "8")
  testthat::expect_equal(dto_Gp, dto$Gp_netMHCpan %>% unique() %>% sort() %>% unique())

  dto_I_pos <- c("0", "1", "2", "4", "5", "6", "7", "8")
  testthat::expect_equal(dto_I_pos, dto$I_pos_netMHC %>% unique() %>% sort() %>% unique())
})
