testthat::test_that("make_BLAST_uuid", {

  if (length(system("which blastp", intern = TRUE)) != 1){
   testthat::skip("Skipping make_BLAST_uuid because ncbiblast+ is not in PATH")
   }

   dt <- data.table::fread("http://get.rech.io/antigen.garnish_test_make_BLAST_uuid.csv")


    dti <- dt %>% data.table::copy %>%
      .[, MHC := "H-2-Kb"] %>%
        .[, nmer_uuid := uuid::UUIDgenerate(), by = "nmer"] %>%
          .[, sample_id := "my_sample"] %>%
            .[, effect_type := "missense_variant"]

    cwd <- getwd()


    Sys.getenv("HOME") %>% setwd

    dto <- make_BLAST_uuid(dti) %>% .[order(nmer)]

    setwd(cwd)

  # run test
       testthat::expect_equal(dto[!is.na(blast_uuid), nmer],
                            c("ENYWRKAY",
                              "ENYWRKSY",
                              "NYWRKAYE",
                              "NYWRKSYE"))

	     testthat::expect_equal(
	        	dto[, pep_type],
	        	c("mutnfs",
							"mutnfs",
							"wt",
							"wt",
							"mutnfs",
							"wt",
							"wt",
							"mutnfs",
						  "wt",
							"wt",
							"wt",
							"mutnfs",
							"wt",
							"mutnfs",
							"wt",
							"wt",
							"mutnfs",
						  "wt",
							"mutnfs",
							"wt"))
    })
