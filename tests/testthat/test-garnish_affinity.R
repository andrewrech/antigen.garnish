
README <- function(){

  testthat::test_that("garnish_affinity strict README example", {
  skip_pred_tools()

# extact copy of README

  # load an example VCF
    dir <- system.file(package = "antigen.garnish") %>%
      file.path(., "extdata/testdata")

    dt <- "antigen.garnish_example.vcf" %>%
    file.path(dir, .) %>%

    # extract variants
      garnish_variants %>%

    # add space separated MHC types
    # see list_mhc() for nomenclature of supported alleles

        .[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")] %>%

    # predict neoantigens
      garnish_affinity

      testthat::expect_true(dt %>% nrow == 483)
      testthat::expect_true(dt[, nmer %>% unique %>% length] == 184)

    })
}

transcripts <- function(){

  testthat::test_that("garnish_affinity from transcripts, diverse MHC", {
  skip_pred_tools()

  # load test data
    dt <- data.table::data.table(
         sample_id = "test",
         ensembl_transcript_id =
         c("ENSMUST00000128119"),
         cDNA_change = c("c.4988C>T"),
         MHC = c("HLA-A*02:01 HLA-DRB1*14:67")) %>%
    # run test
      garnish_affinity(blast = FALSE,
                          save = FALSE)

    testthat::expect_equal(dt[, .N, by = "MHC"]  %>% .[order(MHC)],
    structure(list(MHC = c("HLA-A*02:01", "HLA-DRB1*14:67"),
                    N = c(153L, 30L)),
                    row.names = c(NA, -2L),
                    class = c("data.table", "data.frame"))
                      )
    })
}

excel <- function(){
  testthat::test_that("garnish_affinity from excel file", {

    d <- test_data_dir()

    # load test data
      path <- file.path(d, "antigen.garnish_test_input.xlsx")

    # run test
    dt <- garnish_affinity(path = path,
                              predict = FALSE,
                              blast = FALSE)


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

  d <- test_data_dir()

    # load test data
      dt <- file.path(d, "antigen.garnish_example_jaffa_input.csv") %>%
      data.table::fread

      dt[, MHC := "H-2-Kb"]

    # run test
      dt <- garnish_affinity(dt, blast = FALSE, save = FALSE)

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
                                 blast = FALSE)

    testthat::expect_equal(dto$nmer %>% unique %>% length,
                           109)
      })
  }

  cellular_fraction <- function(){
    testthat::test_that("garnish_affinity with cellular_fraction", {

    skip_pred_tools()

    d <- test_data_dir()

      # load test data
        dt <- file.path(d, "antigen.garnish_test.vcf") %>%

      # run test
        garnish_variants %>%
          .[, MHC := c("HLA-A*02:01")] %>%
          # keep the test small
                      .[1:2] %>%
                      garnish_affinity(save = FALSE)

      testthat::expect_true(dt %>% nrow == 426)
      testthat::expect_true(dt[pep_type %like% "mut",
                              cl_proportion %>%
                              unique %>% signif(digits = 3)] == 0.272)
      testthat::expect_true(dt[pep_type %like% "mut", nmer %>%
                            unique %>% length] == 154)
      testthat::expect_true(dt[!is.na(iedb_score)] %>%
                            nrow == 4)
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
           MHC = "H-2-Kb")

      data.table::data.table(
        id =
        c("ENSMUST00000128119",
        "ENSMUST00000044250",
        "ENSMUST00000018743"),
        test = c(100, 5, 1)) %>%
        data.table::fwrite("antigen.garnish_rna_temp.txt", sep = "\t",
                            quote = FALSE, row.names = FALSE)

      # run test
        dt <- garnish_affinity(dt,
                            counts = "antigen.garnish_rna_temp.txt",
                            min_counts = 10,
                            blast = FALSE,
                            save = FALSE)

      testthat::expect_equal(dt$ensembl_transcript_id %>%
                             unique %>% length, 1)
      testthat::expect_true(all(!dt$ensembl_transcript_id %chin%
                                c("ENSMUST00000018743", "ENSMUST00000044250")))

      })
  }

README()
transcripts()
excel()
jaffa()
peptides()
cellular_fraction()
RNA_test()

list.dirs(path = ".") %>% stringr::str_extract("ag_[a-f0-9]{18}") %>% na.omit %>% lapply(., function(dir){
         unlink(dir, recursive = TRUE, force = TRUE)
        })
