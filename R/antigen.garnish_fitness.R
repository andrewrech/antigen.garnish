## ---- garnish_fitness
#' Run neoantigen fitness model
#'
#' Implement neoantigen fitness model from [Luksza et al. **Nature** 2017](https://www.ncbi.nlm.nih.gov/pubmed/29132144).
#'
#' Source code for this function adapted from example code provided in supplementary data file 7 at this [mirror](https://images.nature.com/original/nature-assets/nature/journal/v551/n7681/extref/nature24473-s9.zip).
#'
#' @param dt A data table created by `garnish_predictions`, to include fitness predictions for non-missense variants, the function must have been called with the argument `blast = TRUE`.
#' @param a Fitness model parameter, see references. `a` represents the horizontal displacement of the binding curve when determining TCR recognition probability of a neoantigen compared to an IEDB near match. Default is 26 as suggested by authors.
#' @param k Fitness model paramater, see references. `k` sets the steepness of the binding curve at `a`.  Yea, I don't get it either. Default is 4.86936 as suggested by authors.
#'
#' @return Returns the input data table with additional columns:
#' * **ResidueChangeClass**: mutant cDNA sequence
#' * **A**: `A` component of fitness model, differential MHC affinity of mutant and closest wt peptide, similar to `BLAST_A`.
#' * **R**: TCR recognition probability, determined by comparison to known epitopes in the IEDB.
#' * **NeoantigenRecognitionPotential**: Product of A and R, the maximum value per sample is the domninant neoepitope.
#'
#' @examples
#'\dontrun{
#'library(antigen.garnish)
#'library(leepR)
#'
#' dt <- "antigen.garnish_example.vcf" %T>%
#' utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#' garnish_variants %>%
#'  .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
#'               "H-2-Kb H-2-IAd",
#'            "HLA-A*01:47 HLA-DRB1*03:08")] %>%
#'              garnish_predictions(blast = TRUE) %>%
#'                garnish_fitness
#'
#'}
#'
#' @import data.table
#' @import dt.inflix
#'
#' @references Luksza, M, Riaz, N, Makarov, V, Balachandran VP, et al. A neoantigen fitness model predicts tumour response to checkpoint blockade immunotherapy **Nature** 2017
#'
#' @export garnish_fitness
#' @md

garnish_fitness <- function(dt,
                            a = 26,
                            k = 4.86936){

          on.exit({
                message("Removing temporary files")
                list.files(pattern = "neoantigens.txt") %>% file.remove
                system("rm -rf fitness_model/")
                              })

  if (dir.exists("Output")) system("sudo rm -rf Output/")

  # check input
  if (class(dt)[1] == "character") dt <- rio::import(dt) %>% data.table::as.data.table

  if (class(dt)[1] == "data.frame") dt %<>% data.table::as.data.table

  if (class(dt)[1] != "data.table") stop("Input must be a data table or file path to garnish_predictions output.")

  if(!(c("var_uuid", "Consensus_scores", "MHC", "nmer", "nmer_uuid", "sample_id") %chin% names(dt)) %>% any) stop("Columns are missing, input must be a data table or file path to garnish_predictions output.")

  if("fus_uuid" %chin% names(dt)) dt[!is.na(fus_uuid) & fus_uuid != "", var_uuid := fus_uuid]

  # construct input table, defer to DAI calculations over BLAST results if it exists

  dta <- dt[pep_type == "wt" & !is.na(dai_uuid), .SD %>% unique,
            .SDcols = c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                        "Consensus_scores", "dai_uuid")] %>%
                data.table::setnames(c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                                       "Consensus_scores", "dai_uuid"),
                                     c("ID", "MUTATION_ID", "Sample", "WT.Peptide", "MT.Allele", "WT.Score",
                                       "link_uuid"))

  dtb <- dt[pep_type != "wt" & !is.na(dai_uuid), .SD %>% unique,
            .SDcols = c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                        "Consensus_scores", "dai_uuid")] %>%
    data.table::setnames(c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                           "Consensus_scores", "dai_uuid"),
                         c("ID", "MUTATION_ID", "Sample", "MT.Peptide", "MT.Allele", "MT.Score",
                           "link_uuid"))

  dti <- merge(dta, dtb, by = c("MUTATION_ID", "Sample", "MT.Allele", "link_uuid"))

  if (nrow(dti) == 0) dti <- data.table::data.table()

  if ("blast_uuid" %chin% names(dt)){

    dta <-  dt[pep_type == "wt" & is.na(dai_uuid) & !is.na(blast_uuid), .SD %>% unique,
                    .SDcols = c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                                "Consensus_scores", "blast_uuid")] %>%
    data.table::setnames(c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                           "Consensus_scores", "blast_uuid"),
                         c("ID", "MUTATION_ID", "Sample", "WT.Peptide", "MT.Allele", "WT.Score",
                           "link_uuid"))

    dtb <-  dt[pep_type != "wt" & is.na(dai_uuid) & !is.na(blast_uuid), .SD %>% unique,
             .SDcols = c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                         "Consensus_scores", "blast_uuid")] %>%
    data.table::setnames(c("nmer_uuid", "var_uuid", "sample_id", "nmer", "MHC",
                           "Consensus_scores", "blast_uuid"),
                         c("ID", "MUTATION_ID", "Sample", "MT.Peptide", "MT.Allele", "MT.Score",
                           "link_uuid"))

  if (nrow(dti) != 0) dti <- data.table::rbindlist(list(dti,
                      merge(dta, dtb,
                    by = c("MUTATION_ID", "Sample", "MT.Allele", "link_uuid"))))

  if (nrow(dti) == 0) dti <-  merge(dta, dtb, by = c("MUTATION_ID", "Sample", "MT.Allele", "link_uuid"))

  }

  # Their code only handles 9 AAs
  dti <- dti[nchar(WT.Peptide) == 9 & nchar(MT.Peptide) == 9]

  if (nrow(dti) == 0) stop("No variants compatible for fitness modeling (mutant and wt 9mers needed).")

  dti <- dti[, ID := link_uuid, by = 1:nrow(dti)]

  # make new IDs compatible with their code
  v <- dti[, ID %>% unique]

  names(v) <- seq_along(v)

    dti[, ID := ID %>% (function(id){
      id <- names(v[which(id == v)])
    }), by = "ID"]

  dti[, MT.Allele := MT.Allele %>% stringr::str_replace_all("-", replacement = "_")]

  dti[, HLA := paste(MT.Allele %>% unique, collapse = " "), by = Sample]

  dti[, chop_score := 1 %>% as.numeric]

  dti[, MUTATION_ID := MUTATION_ID %>% stringr::str_replace_all("-", replacement = "_")]

  dti[, Sample := Sample %>% stringr::str_replace_all("-", replacement = "_")]

  dti[, .SD %>% unique, .SDcols = c("ID", "MUTATION_ID", "Sample", "WT.Peptide", "MT.Peptide",
            "MT.Allele", "WT.Score", "MT.Score", "HLA", "chop_score")] %>%
            data.table::fwrite("neoantigens.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

  # park this for later
  keydt <- dti[, .SD %>% unique, .SDcols = c("MUTATION_ID", "ID.x", "ID.y", "ID", "WT.Peptide", "MT.Peptide")]

  # run blastp
  ## create input

  dtl <- dti %>% split(by = "Sample")

  dir.create("fitness_model")

  lapply(dtl %>% seq_along, function(i){

    dti <- dtl[[i]]

    dti[, label_wt := paste(Sample, "WT", ID, MUTATION_ID, sep = "|")]

    dti[, label_mut := paste(Sample, "MUT", ID, MUTATION_ID, sep = "|")]

    aa <- dti[, WT.Peptide]

    names(aa) <- dti[, label_wt]

    aa_2 <- dti[, MT.Peptide]

    names(aa_2) <- dti[, label_mut]

    aa <- Biostrings::AAStringSet(c(aa, aa_2), use.names = TRUE)

    filename <- paste0("fitness_model/neoantigens_", dti[, Sample %>% unique], ".fasta")

    aa %>% Biostrings::writeXStringSet(filename)

  ## run blastp on fasta

  system(
    paste(
  "blastp -query", filename, "-db /usr/local/bin/iedb.bdb -num_threads", parallel::detectCores(), "-outfmt 5 -evalue 100000000 -gapopen 11 -gapextend 1 >",
  filename %>% stringr::str_replace("\\.fasta", replacement = "_iedb.xml"),
    sep = " ")
  )

  })
  # run pipeline

  if (length(system("which python", intern = TRUE)) != 1) stop("Python is required, please install and rerun.")

  dir.create("Output")

  py_path <- system.file(package = "antigen.garnish") %>% file.path(., "extdata/src/main.py")

  system(
    paste(
  "python", py_path, "neoantigens.txt", "fitness_model", a, k, "Output/neoantigen_fitness_model_output.txt"
    )
  )

  #curate output
  if (!file.exists("Output/neoantigen_fitness_model_output.txt")) stop("Python call did not yield output file.")

  dto <- "Output/neoantigen_fitness_model_output.txt" %>% data.table::fread

  dto[Excluded == FALSE, max(NeoantigenRecognitionPotential), by = Sample] %>% print

  message("Model output table is in Output/")

  dto %>% data.table::setnames(c("Mutation", "Sample"), c("var_uuid", "sample_id"))

  dto %<>% melt(measure.vars = c("WildtypePeptide", "MutantPeptide"), value.name = "nmer")

  # reverse the formatting changes for python scripts

  if ((dt[, sample_id %>% unique] %like% "-") %>% any){

    dt[, sample_id := sample_id %>% stringr::str_replace_all("-", replacement = "_")]

    message("Replacing \"-\" in sample_id with \"_\"...")

  }

  dto[, var_uuid := var_uuid %>% stringr::str_replace_all("_", replacement = "-")]

  dt <- merge(dt, dto[Excluded == FALSE, .SD %>% unique,
          .SDcols = c("var_uuid", "sample_id", "ResidueChangeClass", "A", "R", "NeoantigenRecognitionPotential", "nmer")],
            all.x = TRUE,
              by = c("var_uuid", "sample_id", "nmer"))


  return(dt)

}
