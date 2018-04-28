
README <- function(){

  testthat::test_that("garnish_affinity strict README example", {
  skip_pred_tools()

      # load test data
        dt <- "antigen.garnish_example.vcf" %T>%
        utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

      # run test
        garnish_variants %>%
          .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                       "H-2-Kb H-2-IAd",
                      "HLA-A*01:47 HLA-DRB1*03:08")] %>%
                      garnish_affinity

      testthat::expect_true(dt %>% nrow == 991)
      testthat::expect_true(dt[iedb_score == 1] %>% nrow == 69)
      testthat::expect_true(dt[, nmer %>% unique %>% length] == 574)
    })
}

transcripts <- function(){

  testthat::test_that("garnish_affinity from transcripts, diverse MHC", {
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
      garnish_affinity(blast = FALSE,
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
  testthat::test_that("garnish_affinity from Excel file", {

    # load test data
      path <- "antigen.garnish_test_input.xlsx" %T>%
        utils::download.file("http://get.rech.io/antigen.garnish_test_input.xlsx", .)

    # run test
    dt <- garnish_affinity(path = path,
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
      dt <- garnish_affinity(dt, blast = FALSE, fitness = FALSE)

    testthat::expect_equal(dt %>% nrow, 1071)

     })
}

peptides <- function(){
  testthat::test_that("garnish_affinity assemble from peptides", {

    # load test data
      dt <- data.table::data.table(
              sample_id = "test",
              pep_mut = "AAVMILKWTFRAKINGDEHKRDEF",
              mutant_index = c("7 13 14",  "all", NA),
              MHC = c("H-2-Kb")
                                   )
    # run test data
      dto <- garnish_affinity(dt,
                                 predict = FALSE,
                                 blast = FALSE,
                                 fitness = FALSE)

    testthat::expect_equal(dto$nmer %>% unique %>% length,
                           109)
      })
  }

  cell_fraction <- function(){
    testthat::test_that("garnish_affinity with cell_fraction", {

    skip_pred_tools()

      # load test data
        dt <- "antigen.garnish_test.vcf" %T>%
        utils::download.file("http://get.rech.io/antigen.garnish_test.vcf", .) %>%

      # run test
        garnish_variants %>%
          .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67")] %>%
                      .[1:3] %>%
                      garnish_affinity

      testthat::expect_true(dt %>% nrow == 929)
      testthat::expect_true(dt[pep_type %like% "mut",
                              garnish_score %>%
                              unique %>% signif(digits = 3)] == 1.38)
      testthat::expect_true(dt[, nmer %>%
                            unique %>% length] == 566)
      testthat::expect_true(dt[iedb_score %>% signif(digits = 3) == 1] %>%
                            nrow == 42)
      testthat::expect_true(dt[!is.na(cl_proportion), cl_proportion %>%
                            unique %>% length] == 2)
      })
  }

  RNA_test <- function(){

    testthat::test_that("garnish_affinity from transcripts, with RNA", {
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
           MHC = "H-2-Kb H-2-IAd")

      data.table::data.table(
        id =
        c("ENSMUST00000128119",
        "ENSMUST00000044250",
        "ENSMUST00000018743"),
        test = c(100, 10, 1)) %>%
        data.table::fwrite("antigen.garnish_rna_temp.txt", sep = "\t",
                            quote = FALSE, row.names = FALSE)

      # run test
        dt <- garnish_affinity(dt,
                            counts = "antigen.garnish_rna_temp.txt",
                            blast = FALSE,
                            fitness = FALSE)

      testthat::expect_equal(dt$ensembl_transcript_id %>%
                             unique %>% length, 2)
      testthat::expect_true(all(!dt$ensembl_transcript_id %chin%
                                "ENSMUST00000018743"))

      })
  }

parallel::mclapply(list(README,
									 transcripts,
									 Excel,
									 jaffa,
									 peptides,
									 cell_fraction,
									 RNA_test),
									 test_runner)
