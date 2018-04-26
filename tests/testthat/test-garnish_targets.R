testthat::test_that("garnish_antigens",
{
skip_pred_tools()

    # load test data
      dt <- "http://get.rech.io/antigen.garnish_PureCN_example_output.txt" %>%
      data.table::fread %>%
      garnish_antigens

    # run test

    testthat::expect_equal(dt, structure(
     list(sample_id = "antigen.garnish_test_pureCN.vcf",
					nmer = "LICPYMEPI",
					MHC = "HLA-A*02:01",
					external_gene_name = "MTOR",
					protein_change = "p.Arg769Cys",
					cDNA_change = "c.2305C>T",
					Consensus_scores = 132.615664199976,
					fitness_score = 12.5057569549273,
					iedb_score = 0.999999999973312,
					min_DAI = 12.505756955261,
					clone_id = 2L,
					cl_proportion = 0.277452002724779),
					.Names = c("sample_id",
					"nmer",
					"MHC",
					"external_gene_name",
					"protein_change",
					"cDNA_change",
					"Consensus_scores",
					"fitness_score",
					"iedb_score",
					"min_DAI",
					"clone_id",
					"cl_proportion"),
					class = c("data.table",
					"data.frame"), row.names = c(NA, -1L)))

})
