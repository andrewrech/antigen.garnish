on.exit({

  lapply(c("antigen.garnish_example.vcf",
  "antigen.garnish_example_mutect2.vcf",
  "antigen.garnish_example_strelka.vcf",
  "antigen.garnish_GRCh38_pep.RDS",
  "antigen.garnish_GRCm38_pep.RDS",
  "Hu_nmer_fasta.fa",
  "hublastpout.csv",
  "iedb_query.fa",
  "iedbout.csv",
  "Ms_nmer_fasta.fa",
  "msblastpout.csv",
  "neoantigens_fitness_model_output.txt"), function(i){

    if (file.exists(i))
      file.remove(i)
  })
})

if(requireNamespace("testthat", quietly = TRUE)){
    library(testthat)
    library(antigen.garnish)
    # debug(test_file)
    # debug(testthat:::test_dir)
    test_check(package = "antigen.garnish",
               wrap = FALSE,
               stop_on_failure = FALSE,
               filter = "testthat_")
}

