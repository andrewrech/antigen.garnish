
README <- function(){

  testthat::test_that("garnish_predictions strict README example", {
  skip_pred_tools()

      # load test data
        dt <- "antigen.garnish_example.vcf" %T>%
        utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

      # run test
        garnish_variants %>%
          .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                       "H-2-Kb H-2-IAd",
                      "HLA-A*01:47 HLA-DRB1*03:08")] %>%
                      garnish_predictions

      testthat::expect_true(dt %>% nrow == 911)
      testthat::expect_true(dt[iedb_score == 1] %>% nrow == 69)
      testthat::expect_true(dt[, nmer %>% unique %>% length] == 566)
    })
}

transcripts <- function(){

  testthat::test_that("garnish_predictions from transcripts, diverse MHC", {
  skip_pred_tools()

  # load test data
    dt <- data.table::data.table(
         sample_id = "test",
         ensembl_transcript_id =
         c("ENSMUST00000128119",
           "ENSMUST00000044250",
           "ENSMUST00000018743"),
         cDNA_change = c("c.4988C>T",
                         "c.1114T>G",
                         "c.718T>A"),
         MHC = c("HLA-A*02:01 HLA-DRB1*14:67",
                 "H-2-Kb H-2-IAd",
                 "HLA-C*14:02")) %>%
    # run test
      garnish_predictions(blast = FALSE,
                          fitness = FALSE)

    testthat::expect_equal(dt$cDNA_change %>% unique %>% length, 3)
    testthat::expect_equal(dt$MHC %>% unique %>% sort,
        c("H-2-IAd",
          "H-2-Kb",
          "HLA-A*02:01",
          "HLA-C*14:02",
          "HLA-DRB1*14:67"))
    })
}

Excel <- function(){
  testthat::test_that("garnish_predictions from Excel file", {

    # load test data
      path <- "antigen.garnish_test_input.xlsx" %T>%
        utils::download.file("http://get.rech.io/antigen.garnish_test_input.xlsx", .)

    # run test
    dt <- garnish_predictions(path = path,
                              predict = FALSE,
                              blast = FALSE,
                              fitness = FALSE)


    testthat::expect_equal(dt %>% nrow, 551)
    testthat::expect_equal(
      dt[, nmer %>%
           nchar %>%
           unique] %>% sort,
      8:15)

  })
}

jaffa <- function(){
  testthat::test_that("test predict from jaffa input", {

  skip_pred_tools()

    # load test data
      dt <- data.table::fread("http://get.rech.io/antigen.garnish_example_input_jaffa.csv")

      dt[, MHC := "H-2-Kb"]

    # run test
      dt <- garnish_predictions(dt, blast = FALSE, fitness = FALSE)

    testthat::expect_equal(dt %>% nrow, 1071)

     })
}

peptides <- function(){
  testthat::test_that("garnish_predictions assemble from peptides", {

    # load test data
      dt <- data.table::data.table(
              sample_id = "test",
              pep_mut = "AAVMILKWTFRAKINGDEHKRDEF",
              mutant_index = c("7 13 14",  "all", NA),
              MHC = c("H-2-Kb")
                                   )
    # run test data
      dto <- garnish_predictions(dt,
                                 predict = FALSE,
                                 blast = FALSE,
                                 fitness = FALSE)

    testthat::expect_equal(dto$nmer %>% unique %>% length,
                           109)
      })
  }

  CELLFRACTION <- function(){
    testthat::test_that("garnish_predictions with CELLFRACTION", {
    skip_pred_tools()

        # load test data
          dt <- "antigen.garnish_test_pureCN.vcf" %T>%
          utils::download.file("http://get.rech.io/antigen.garnish_test_pureCN.vcf", .) %>%

        # run test
          garnish_variants %>%
            .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67")] %>%
                        .[1:3] %>%
                        garnish_predictions

        testthat::expect_true(dt %>% nrow == 929)
        testthat::expect_true(dt[pep_type %like% "mut",
                                garnish_score %>% unique %>% signif(digits = 3)] == 1.38)
        testthat::expect_true(dt[, nmer %>% unique %>% length] == 566)
        testthat::expect_true(dt[iedb_score %>% signif(digits = 3) == 1] %>% nrow == 42)
        testthat::expect_true(dt[!is.na(clone_prop), clone_prop %>% unique %>% length] == 2)
        
      })
  }

parallel::mclapply(list(README, transcripts, Excel, jaffa, peptides, CELLFRACTION), test_runner)
