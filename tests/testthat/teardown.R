
# tear down antigen.garnish tests

  lapply(
         c(
      "antigen.garnish_example.vcf",
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
      "antigen.garnish_test_input.xlsx",
     "9_neoantigens_Lukza_model_output.txt",
     "antigen.garnish_test_pureCN.vcf",
     "antigen.garnish_pureCN_example_output.txt",
      "antigen.garnish_rna_temp.txt",
      "antigen.garnish_prop_AF.csv",
      "antigen.garnish_prop_CF.csv"), function(i){

        if (file.exists(i))
          file.remove(i)
    })

list.dirs(path = ".") %>% stringr::str_extract("[a-f0-9]{8}-[a-f0-9]{4}-4[a-f0-9]{3}-[89ab][a-f0-9]{3}-[a-f0-9]{12}") %>% na.omit %>% lapply(., function(dir){
        unlink(dir, recursive = TRUE, force = TRUE)
        })
