library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)
library(Biostrings)

testthat::test_that("garnish_predictions", {
  
  if (!check_pred_tools() %>% unlist %>% all) {
    testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
  }
  
  # copy of README
  
  # download an example VCF
  dt <- "antigen.garnish_example.vcf" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
    
    # extract variants
    garnish_variants %>%
    
    # add MHC types
    .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                 "HLA-A*02:01 HLA-DRB1*14:67",
                 "HLA-A*03:01 HLA-DRB1*03:01")] %>%
    
    # predict neoepitopes
    garnish_predictions
  
  # test if nmers are contained in peptides
  dt[pep_type != "wt", test := str_extract(pattern = nmer, string = pep_mut)]
  dt[pep_type == "wt", test := str_extract(pattern = nmer, string = pep_wt)]
  dt[, validate := test == nmer, by = 1:nrow(dt)]
  testthat::expect_equal(dt$validate %>% all(), TRUE)
  
  # test that DAI peptides are only 1 AA apart
  DAIdt <- dt[!is.na(dt$DAI)][pep_type != "wt"]
  DAIdt[, n.aa.mismatch := lapply(1:nrow(DAIdt), function(i){
    x <- pep_mut[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    y <- pep_wt[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    return(sum(x != y))
  }) %>% unlist()] 
  testthat::expect_equal(all(DAIdt$n.aa.mismatch == 1), TRUE)
  
  # match mutant_loc to SnpEff protein change call for all mutants
  dt[, prot.change.from.se := str_extract_all(string = protein_change, pattern =  "[0-9]+") %>% unlist()]
  testthat::expect_equal(dt[, mutant_loc == prot.change.from.se] %>% all(), TRUE)
  
  # test that pep_mut matches SnpEff protein change call for missense
  
  DAIdt[, new_aa_pep := lapply(1:nrow(DAIdt), function(i){
    x <- pep_mut[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    y <- pep_wt[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    nAA <- which(x != y)
    return(x[nAA])
  }) %>% unlist()] 
  DAIdt[, old_aa_pep := lapply(1:nrow(DAIdt), function(i){
    x <- pep_mut[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    y <- pep_wt[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    nAA <- which(x != y)
    return(y[nAA])
  }) %>% unlist()]
  
  convert.aa <- function(x){
    x <- Biostrings::AMINO_ACID_CODE[which(names(Biostrings::AMINO_ACID_CODE) == x)]
    unlist(x)
  }
  
  DAIdt[, old_aa_pep := convert.aa(old_aa_pep), by = 1:nrow(DAIdt)][, new_aa_pep := convert.aa(new_aa_pep), by = 1:nrow(DAIdt)]
  
  DAIdt[, old_aa_se := str_extract(string = protein_change, pattern = "(?<=(p\\.))[A-Z][a-z]{2}(?=[0-9])") %>% unlist()]
  DAIdt[, new_aa_se := str_extract(string = protein_change, pattern = "(?<=[0-9])[A-Z][a-z]{2}$") %>% unlist()]
  
  testthat::expect_equal(DAIdt[, old_aa_pep == old_aa_se] %>% all(), TRUE)
  testthat::expect_equal(DAIdt[, new_aa_pep == new_aa_se] %>% all(), TRUE)
  
  # check wt and mut peptides are registered correctly
  
  ## getting missense AA change
  
  DAIdt[, new_aa_pep := lapply(1:nrow(DAIdt), function(i){
    x <- pep_mut[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    y <- pep_wt[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    nAA <- which(x != y)
    return(x[nAA])
  }) %>% unlist()] 
  DAIdt[, old_aa_pep := lapply(1:nrow(DAIdt), function(i){
    x <- pep_mut[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    y <- pep_wt[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    nAA <- which(x != y)
    return(y[nAA])
  }) %>% unlist()]
  
  ## deriving wt from mutant peptide
  
  DAIdt[, pep_wt2 := lapply(1:nrow(DAIdt), function(i){
    x <- pep_mut[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
    if (x[mutant_loc[i]] == new_aa_pep[i]){
      x[mutant_loc[i]] <- old_aa_pep[i]
    } else{return(NA)}
    return(paste(x, collapse = ""))}) %>% unlist()] 
 
 testthat::expect_equal(DAIdt[, pep_mut == pep_wt2] %>% all(), TRUE)
 
 ## deriving mutant from wt peptide
 
 DAIdt[, pep_mut2 := lapply(1:nrow(DAIdt), function(i){
   x <- pep_wt[i] %>% strsplit(., split = "", fixed = TRUE) %>% unlist()
   if (x[mutant_loc[i]] == old_aa_pep[i]){
     x[mutant_loc[i]] <- new_aa_pep[i]
   } else{return(NA)}
   return(paste(x, collapse = ""))}) %>% unlist()] 

testthat::expect_equal(DAIdt[, pep_mut == pep_mut2] %>% all(), TRUE)

# Compare actual to expected number of peptides generated

ndt <- DAIdt[, c("nmer_uuid", "var_uuid")] %>% unique()
if (any(ndt[, .N, by = var_uuid]$N < 77)){
ndt[, less := .N < 77, by = var_uuid]
vect <- ndt[less == TRUE]$var_uuid
test_dt <- DAIdt[var_uuid %chin% vect]
test_dt[, nchar := nchar(pep_mut)]
test_dt[, min := min(c(mutant_loc, (nchar - mutant_loc))), by = 1:nrow(test_dt)]
testthat::expect_equal(all(test_dt$min < 14), TRUE) 
}
  if (file.exists("antigen.garnish_example.vcf"))
    file.remove("antigen.garnish_example.vcf")
})