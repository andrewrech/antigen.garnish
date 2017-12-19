testthat::test_that("write_mhcnuggets_nmers", {

  list.files(pattern = "mhcnuggets_input_.*csv") %>% file.remove
  on.exit(list.files(pattern = "mhcnuggets_input_.*csv") %>% file.remove)

  # load test data
    alleles <- data.table::rbindlist(
                           list(
    system.file("extdata",
        "netMHC_alleles.txt", package = "antigen.garnish") %>%
                        data.table::fread(header = FALSE, sep = "\t") %>%
                        data.table::setnames("V1", "allele") %>%
                        .[, type := "netMHC"],
    system.file("extdata",
        "netMHCpan_alleles.txt", package = "antigen.garnish") %>%
                        data.table::fread(header = FALSE, sep = "\t") %>%
                        data.table::setnames("V1", "allele") %>%
                        .[, type := "netMHCpan"],
    system.file("extdata",
        "mhcflurry_alleles.txt", package = "antigen.garnish") %>%
                        data.table::fread(header = FALSE, sep = "\t") %>%
                        data.table::setnames("V1", "allele") %>%
                        .[, type := "mhcflurry"],
    system.file("extdata",
        "netMHCII_alleles.txt", package = "antigen.garnish") %>%
                        data.table::fread(header = FALSE, sep = "\t") %>%
                        data.table::setnames("V1", "allele") %>%
                        .[, type := "netMHCII"],
    system.file("extdata",
        "netMHCIIpan_alleles.txt", package = "antigen.garnish") %>%
                        data.table::fread(header = FALSE, sep = "\t") %>%
                        data.table::setnames("V1", "allele") %>%
                        .[, type := "netMHCIIpan"],
    system.file("extdata",
         "mhcnuggets_gru_alleles.txt", package = "antigen.garnish") %>%
                        data.table::fread(header = FALSE, sep = "\t") %>%
                        data.table::setnames("V1", "allele") %>%
                        .[, type := "mhcnuggets_gru"],
    system.file("extdata",
        "mhcnuggets_lstm_alleles.txt", package = "antigen.garnish") %>%
                        data.table::fread(header = FALSE, sep = "\t") %>%
                        data.table::setnames("V1", "allele") %>%
                        .[, type := "mhcnuggets_lstm"]
                        ))

    dt <- data.table::data.table(
       mhcnuggets = c("HLA-A0201",
                      "HLA-A0201",
                      "HLA-A0201",
                      "HLA-A0201",
                      "HLA-A0201",
                      "HLA-A0201"),
       nmer = c("AQSGTPPT",
                "AQSGTPPTG",
                "AQSGTPPTGL",
                "AQSGTPPT",
                "AQSGTPPTG",
                "AQSGTPPTGL"),
       nmer_l = c(8L, 9L, 10L, 8L, 9L, 10L))

    # run test
      write_mhcnuggets_nmers(dt, alleles)

      out <- list.files(pattern = "mhcnuggets_input_.*csv")

      testthat::expect_equal(
        out %>% length, 2)

      testthat::expect_equal(

out %>% lapply(., function(x){
          x %>% fread(header = FALSE)
        }) %>% unlist %>% length, 6)

    })
