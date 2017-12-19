testthat::test_that("extract_cDNA", {

  # load test data
    dt <- data.table::data.table(cDNA_change =
        c("c.4988C>T",
          "c.1114T>G",
          "c.718T>A",
          "c.718_722delTTGC",
          "c.315_316insTAATGATACGGC",
          "c.3088_3089insT")) %>%
  # run test
    extract_cDNA

  testthat::expect_equal(dt, structure(list(cDNA_change = c("c.4988C>T", "c.1114T>G", "c.718T>A",
    "c.718_722delTTGC", "c.315_316insTAATGATACGGC", "c.3088_3089insT"),
    cDNA_locs = c(4988L, 1114L, 718L, 718L, 315L, 3088L),
    cDNA_locl = c(4988L, 1114L, 718L, 722L, 316L, 3089L),
    cDNA_type = c(">", ">", ">", "del", "ins", "ins"),
    cDNA_seq = c("T", "G", "A", "TTGC", "TAATGATACGGC", "T")),
    .Names = c("cDNA_change", "cDNA_locs", "cDNA_locl", "cDNA_type", "cDNA_seq"),
    row.names = c(NA, -6L),
    class = c("data.table",
    "data.frame"))
)
    })