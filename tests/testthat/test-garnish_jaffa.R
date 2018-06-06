testthat::test_that("garnish_jaffa", {

    d <- test_data_dir()

    list.files(pattern = "antigen.garnish_jaffa") %>% file.remove
    on.exit(list.files(pattern = "antigen.garnish_jaffa") %>% file.remove)

  # load test data
    path <- file.path(d, "antigen.garnish_jaffa_results.csv")
    fasta_path <- file.path(d, "antigen.garnish_jaffa_results.fasta")

  # run test
   dt <- garnish_jaffa(path = path, fasta_path = fasta_path)

  testthat::expect_equal(dt %>% nrow, 15)
  testthat::expect_equal(dt %>% length, 11)
    })
