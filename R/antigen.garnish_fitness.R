## ---- garnish_fitness
#' Generate neoepitope fitness model
#'
#' Implements the neoepitopes fitness model of [Luksza et al. *Nature* 2017](https://www.ncbi.nlm.nih.gov/pubmed/29132144).
#'
#' @param dt Data table output from `garnish_predictions`.
#' @param a Fitness model parameter. Binding curve horizontal displacement used to determine TCR recognition probability of a peptide compared to an IEDB near match.
#' @param k Fitness model parameter. Steepness of the binding curve at `a`
#'
#' @return A data table with additional columns:
#' * **ResidueChangeClass**: mutant cDNA sequence
#' * **A**: `A` component of fitness model, differential MHC affinity of mutant and closest wt peptide, similar to `BLAST_A`.
#' * **R**: TCR recognition probability, determined by comparison to known epitopes in the IEDB.
#' * **NeoantigenRecognitionPotential**: Product of A and R, the maximum value per sample is the dominant neoepitope.
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
#' @references Luksza, M, Riaz, N, Makarov, V, Balachandran VP, et al. A neoepitopes fitness model predicts tumour response to checkpoint blockade immunotherapy **Nature** 2017
#'
#' @export garnish_fitness
#' @md

garnish_fitness <- function(dt,
                            a = 26,
                            k = 4.86936){

    on.exit({
          message("Removing temporary files")
          list.files(pattern = "neoepitopes.txt") %>% file.remove
          unlink("fitness_model", recursive = TRUE, force = TRUE)
                              })

  # check input

    if ("data.table" %chin% class(dt))
      stop("Input must be a data table.")

    if(!(c("var_uuid",
           "Consensus_scores",
           "MHC",
           "nmer",
           "nmer_uuid",
           "sample_id") %chin% names(dt)) %>% any)
      stop("
Input data table must include
the following columns:

  var_uuid
  Consensus_scores
  MHC
  nmer
  nmer_uuid
  sample_id
      ")

    if("fus_uuid" %chin% names(dt))
      dt[!is.na(fus_uuid) & fus_uuid != "", var_uuid := fus_uuid]

  # construct input table, defer to DAI calculations over BLAST results if it exists

            dta <- dt[pep_type == "wt" & !is.na(dai_uuid), .SD %>% unique,
                      .SDcols = c("nmer_uuid",
                                  "var_uuid",
                                  "sample_id",
                                  "nmer",
                                  "MHC",
                                  "Consensus_scores",
                                  "dai_uuid")] %>%
            data.table::setnames(c("nmer_uuid",
                                   "var_uuid",
                                   "sample_id",
                                   "nmer",
                                   "MHC",
                                   "Consensus_scores",
                                   "dai_uuid"),
                                 c("ID",
                                   "MUTATION_ID",
                                   "Sample",
                                   "WT.Peptide",
                                   "MT.Allele",
                                   "WT.Score",
                                   "link_uuid"))

            dtb <- dt[pep_type != "wt" & !is.na(dai_uuid), .SD %>% unique,
                      .SDcols = c("nmer_uuid",
                                  "var_uuid",
                                  "sample_id",
                                  "nmer",
                                  "MHC",
                                  "Consensus_scores",
                                  "dai_uuid")] %>%
              data.table::setnames(
                                c("nmer_uuid",
                                  "var_uuid",
                                  "sample_id",
                                  "nmer",
                                  "MHC",
                                  "Consensus_scores",
                                  "dai_uuid"),
                                c("ID",
                                  "MUTATION_ID",
                                  "Sample",
                                  "MT.Peptide",
                                  "MT.Allele",
                                  "MT.Score",
                                  "link_uuid"))

    dti <- data.table::copy(merge(dta, dtb, by = c("MUTATION_ID",
                                      "Sample",
                                      "MT.Allele",
                                      "link_uuid")))

    if (nrow(dti) == 0) dti <- data.table::data.table()

    if ("blast_uuid" %chin% names(dt)){

      dta <-  dt[pep_type == "wt" & is.na(dai_uuid) & !is.na(blast_uuid), .SD %>% unique,
                      .SDcols = c("nmer_uuid",
                                  "var_uuid",
                                  "sample_id",
                                  "nmer",
                                  "MHC",
                                  "Consensus_scores",
                                  "blast_uuid")] %>%
           data.table::setnames(c("nmer_uuid",
                                  "var_uuid",
                                  "sample_id",
                                  "nmer",
                                  "MHC",
                                  "Consensus_scores",
                                  "blast_uuid"),
                                c("ID",
                                  "MUTATION_ID",
                                  "Sample",
                                  "WT.Peptide",
                                  "MT.Allele",
                                  "WT.Score",
                                  "link_uuid"))

             dtb <-  dt[pep_type != "wt" & is.na(dai_uuid) & !is.na(blast_uuid), .SD %>% unique,
                      .SDcols = c("nmer_uuid",
                                  "var_uuid",
                                  "sample_id",
                                  "nmer",
                                  "MHC",
                                  "Consensus_scores",
                                  "blast_uuid")] %>%
             data.table::setnames(
                                c("nmer_uuid",
                                  "var_uuid",
                                  "sample_id",
                                  "nmer",
                                  "MHC",
                                  "Consensus_scores",
                                  "blast_uuid"),
                                c("ID",
                                  "MUTATION_ID",
                                  "Sample",
                                  "MT.Peptide",
                                  "MT.Allele",
                                  "MT.Score",
                                  "link_uuid"))

    if (nrow(dti) != 0)
      dti <- data.table::copy(
              data.table::rbindlist(list(dti,
                merge(dta, dtb,
                      by = c("MUTATION_ID",
                             "Sample",
                             "MT.Allele",
                             "link_uuid")))))

    if (nrow(dti) == 0)
        dti <- merge(dta, dtb,
                     by = c("MUTATION_ID",
                            "Sample",
                            "MT.Allele",
                            "link_uuid"))
}

  # use only unique 9mers

    dti <- dti[nchar(WT.Peptide) == 9 & nchar(MT.Peptide) == 9]

    if (nrow(dti) == 0) {
      warning("No 9mers compatible for fitness modeling.")
      return(dt)
    }

    dti <- dti[, ID := link_uuid, by = 1:nrow(dti)]

    # make compatible IDs
      v <- dti[, ID %>% unique]
      names(v) <- seq_along(v)

        dti[, ID := ID %>% (function(id){
          id <- names(v[which(id == v)])
        }), by = "ID"]

  for (col in c("MT.Allele",
                "MUTATION_ID",
                "Sample"))
  data.table::set(dti, j = col,
                  value = dti[[eval(col)]] %>%
                  stringr::str_replace_all(stringr::fixed("-"), replacement = "_"))

  dti[, HLA := paste(MT.Allele %>%
                      unique, collapse = " "), by = Sample]
  dti[, chop_score := 1]
  dti[, .SD %>% unique,
    .SDcols = c("ID",
                "MUTATION_ID",
                "Sample",
                "WT.Peptide",
                "MT.Peptide",
                "MT.Allele",
                "WT.Score",
                "MT.Score",
                "HLA",
                "chop_score")] %>%
            data.table::fwrite("neoepitopes.txt",
                                sep = "\t",
                                quote = FALSE,
                                col.names = TRUE,
                                row.names = FALSE)

  # run blastp
    # create input

      dtl <- dti %>% split(by = "Sample")

      dir.create("fitness_model", showWarnings = FALSE)

      lapply(dtl %>% seq_along, function(i){
        dti <- data.table::copy(dtl[[i]])
        dti[, label_wt := paste(Sample, "WT", ID, MUTATION_ID, sep = "|")]
        dti[, label_mut := paste(Sample, "MUT", ID, MUTATION_ID, sep = "|")]
        aa <- dti[, WT.Peptide]
        names(aa) <- dti[, label_wt]
        aa_2 <- dti[, MT.Peptide]
        names(aa_2) <- dti[, label_mut]
        aa <- Biostrings::AAStringSet(c(aa, aa_2), use.names = TRUE)
        filename <- paste0("fitness_model/neoepitopes_", dti[, Sample %>% unique], ".fasta")
        aa %>% Biostrings::writeXStringSet(filename)

    # run blastp on fasta

      system(
        paste(
      "blastp -query", filename, "-db ~/antigen.garnish/iedb.bdb -num_threads", parallel::detectCores(), "-outfmt 5 -evalue 100000000 -gapopen 11 -gapextend 1 >",
      filename %>% stringr::str_replace("\\.fasta", replacement = "_iedb.xml"),
        sep = " ")
      )

    })

  # run pipeline

  py_path <- system.file(package = "antigen.garnish") %>% file.path(., "extdata/src/main.py")

  system(
    paste("python",
          py_path,
          "neoepitopes.txt",
          "fitness_model",
          a,
          k,
          "neoepitopes_fitness_model_output.txt"
    )
  )

  # curate output
  if (!file.exists("neoepitopes_fitness_model_output.txt")) {
    warning("garnish_fitness did not return output.")
    return(dt)
  }

  dto <- "neoepitopes_fitness_model_output.txt" %>% data.table::fread

  dto %>% data.table::setnames(c("Mutation", "Sample"),
                               c("var_uuid", "sample_id"))
  dto %<>% melt(measure.vars = c("WildtypePeptide", "MutantPeptide"),
                value.name = "nmer")

  # remove existing columns
    dt[, .SD := NULL,
         .SDcols = c(
                     "ResidueChangeClass",
                     "A",
                     "R",
                     "NeoantigenRecognitionPotential"
                    )]

  dto[, var_uuid := var_uuid %>%
        stringr::str_replace_all(stingr::fixed("_"), replacement = "-")]

  dt <- merge(dt, dto[Excluded == FALSE, .SD %>% unique,
          .SDcols = c("var_uuid",
                      "sample_id",
                      "ResidueChangeClass",
                      "A",
                      "R",
                      "NeoantigenRecognitionPotential",
                      "nmer")],
                      all.x = TRUE,
                      by = c("var_uuid",
                             "sample_id",
                             "nmer"))

  return(dt)

}
