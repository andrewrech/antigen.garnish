testthat::test_that("get_vcf_snpeff_dt", {

  # load test data

		dt <- "antigen.garnish_test.vcf" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_test.vcf", .) %>%
    	vcfR::read.vcfR(vcfo, verbose = TRUE) %>%
      get_vcf_sample_dt

  # run test

    dt[, .SD, .SDcols = c("effect_type", "ensembl_transcript_id",
    "ensembl_gene_id", "protein_change", "cDNA_change",
    "protein_coding")] %>%

    testthat::expect_equal(dt %>% names,
                           c("FORMAT",
														"normal_sample",
														"tumor_sample",
														"normal_sample_GT",
														"normal_sample_BQ",
														"normal_sample_DP",
														"normal_sample_FA",
														"normal_sample_SS",
														"tumor_sample_GT",
														"tumor_sample_BQ",
														"tumor_sample_DP",
														"tumor_sample_FA",
														"tumor_sample_SS",
														"normal_sample_AD_ref",
														"normal_sample_AD_alt",
														"tumor_sample_AD_ref",
														"tumor_sample_AD_alt"))

    testthat::expect_equal(dt %>% nrow, 61)
    testthat::expect_equal(dt %>% length, 17)

    })