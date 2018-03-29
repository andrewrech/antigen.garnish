testthat::test_that("garnish_variants intersection", {

  # load test data
    "antigen.garnish_example_mutect2.vcf" %>%
      utils::download.file(paste0("http://get.rech.io/", .), .)

    "antigen.garnish_example_strelka.vcf" %>%
      utils::download.file(paste0("http://get.rech.io/", .), .)

  # run test
    dto_i <- garnish_variants(c("antigen.garnish_example_mutect2.vcf",
                     "antigen.garnish_example_strelka.vcf"))

  testthat::expect_equal(dto_i %>% length, 23)
  testthat::expect_equal(dto_i %>% nrow, 32)

  # run test
  dto_ni <- garnish_variants(c("antigen.garnish_example_mutect2.vcf",
                     "antigen.garnish_example_strelka.vcf"),
            intersect = FALSE)

  testthat::expect_equal(dto_ni %>% length, 25)
  testthat::expect_equal(dto_ni %>% nrow, 222)

    })