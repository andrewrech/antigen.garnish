
run_path_error <- function() {
  testthat::test_that("antigen.garnish_run path error", {
    err <- try(
      {
        Sys.setenv(AG_DATA_DIR = "fake_directory")

        dt <- data.table::data.table(
          sample_id = "test",
          transcript_id =
            c("ENSMUST00000128119"),
          cDNA_change = c("c.4988C>T"),
          MHC = c("HLA-A*02:01 HLA-DRB1*14:67")
        ) %>%

          # run test
          garnish_affinity()
      },
      silent = TRUE
    )

    testthat::expect_true("try-error" == (err %>% class()))

    Sys.setenv(AG_DATA_DIR = "")
  })
}

error_function <- function() {
  testthat::test_that("ag error report", {
    err <- try(
      {
        .ag_data_err()
      },
      silent = TRUE
    )

    testthat::expect_true("try-error" == (err %>% class()))
  })
}

run_path_error()
error_function()
