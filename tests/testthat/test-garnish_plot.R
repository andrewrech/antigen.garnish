testthat::test_that("garnish_plot", {

  if (identical(Sys.getenv("TESTTHAT"), "true"))
    on.exit(list.files(pattern = "(antigen.garnish.*\\.pdf)|(antigen.garnish_example_output)") %>% file.remove)

    d <- test_data_dir()

    dt <- data.table::fread(file.path(d, "antigen.garnish_PureCN_example_output.txt"))

    garnish_plot(dt)

  # run test
    testthat::expect_equal(list.files(pattern = "antigen.garnish.*summary.*\\.pdf") %>% length, 1)
    testthat::expect_equal(list.files(pattern = "antigen.garnish.*score.*\\.pdf") %>% length, 1)

})
