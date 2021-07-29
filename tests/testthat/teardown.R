
# tear down tests

lapply(
  c(
    "antigen.garnish_example.vcf",
    "antigen.garnish_GRCh38_pep.RDS",
    "antigen.garnish_GRCm38_pep.RDS",
    "Hu_ag_nmer_fasta.fa",
    "hublastpout.csv",
    "iedb_query.fa",
    "iedbout.csv",
    "Ms_ag_nmer_fasta.fa",
    "msblastpout.csv",
    "antigen.garnish_test_input.xlsx",
    "antigen.garnish_collate_example.txt",
    "antigen.garnish_example_output.txt",
    "antigen.garnish_rna_temp.txt"
  ), function(i) {
    if (file.exists(i)) {
      file.remove(i)
    }
  }
)

list.files(pattern = "^Luksza_model|^ag_output|^a\\.g_ex.*csv$") %>% unlink(force = TRUE)

list.dirs(path = ".") %>%
  stringr::str_extract(".*[a-f0-9]{18}") %>%
  na.omit() %>%
  lapply(., function(dir) {
    unlink(dir, recursive = TRUE, force = TRUE)
  })
