testthat::test_that("get_vcf_snpeff_dt", {

  # load test data

		dt <- "antigen.garnish_test.vcf" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_test.vcf", .) %>%
    	vcfR::read.vcfR(verbose = TRUE) %>%
      get_vcf_info_dt

    testthat::expect_equal(dt %>% names,
   												c(
												   "CHROM",
												   "POS",
												   "ID",
												   "REF",
												   "ALT",
												   "QUAL",
												   "FILTER",
												   "INFO",
												   "AC",
												   "AF",
												   "AN",
												   "DB",
												   "SOMATIC",
												   "VT",
												   "set",
												   "OnTarget",
												   "SM0",
												   "SM1",
												   "SM2",
												   "SM3",
												   "SM4",
												   "SM5",
												   "SM6",
												   "SM7",
												   "GM0",
												   "GM1",
												   "GM2",
												   "GM3",
												   "GM4",
												   "GM5",
												   "GM6",
												   "GM7",
												   "GCONTHIGH",
												   "GCONTLOW",
												   "GHOMOZYGOUS",
												   "ML.SOMATIC",
												   "ML.M",
												   "ML.C",
												   "ML.M.SEGMENT",
												   "ML.AR",
												   "CF",
												   "PS",
												   "LR",
												   "ANN",
												   "Cosmic.CNT",
												   "LOF",
												   "ML.LOH"
   											  	))

    testthat::expect_equal(dt %>% nrow, 61)
    testthat::expect_equal(dt %>% length, 47)

    })