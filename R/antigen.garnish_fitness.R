## ---- garnish_fitness
#' Generate neoepitope fitness model
#'
#' Implements the neoepitope fitness model of [Luksza et al. *Nature* 2017](https://www.ncbi.nlm.nih.gov/pubmed/29132144).
#'
#' @param dt Data table output from `garnish_predictions`.
#' @param a Fitness model parameter. Binding curve horizontal displacement used to determine TCR recognition probability of a peptide compared to an IEDB near match.
#' @param k Fitness model parameter. Steepness of the binding curve at `a`.
#'
#' @return A data table with added fitness model parameter columns:
#' * **ResidueChangeClass**: Amino acid change class, eg hydrophobic to non-hydrophobic.
#' * **R**: TCR recognition probability, determined by comparison to known epitopes in the IEDB.
#' * **NeoantigenRecognitionPotential**: Neoantigen Recognition Potential calculated by the model. Original validation was limited to 9mers, class I, and human IEDB epitopes. The model will be applied to all peptides, class I and II.
#'
#' @export garnish_fitness
#' @md

garnish_fitness <- function(dt,
                            a = 26,
                            k = 4.86936){

    on.exit({
          message("Removing temporary files")
          try(
          list.files(pattern = "neoepitopes.txt") %>% file.remove
          )
          try(
          list.files(pattern = "fitness_model_output.txt") %>% file.remove
          )
          try(
          unlink(list.files(pattern = "fitness_model_[0-9]+"), recursive = TRUE, force = TRUE)
          )
                              })

  # check input

    if (!"data.table" %chin% class(dt))
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

            dta <- dt[pep_type == "wt" & !is.na(dai_uuid) & dai_uuid != "", .SD %>% unique,
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

            dtb <- dt[pep_type != "wt" & !is.na(dai_uuid) & dai_uuid != "", .SD %>% unique,
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

      dta <-  dt[pep_type == "wt" & (is.na(dai_uuid) | dai_uuid == "") &
                  !is.na(blast_uuid) & blast_uuid != "", .SD %>% unique,
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

             dtb <-  dt[pep_type != "wt" & (is.na(dai_uuid) | dai_uuid == "") &
                         !is.na(blast_uuid) & blast_uuid != "", .SD %>% unique,
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

  # loop over nmer lengths and species
  dti[MT.Allele %like% "H-2", spc := "Ms"] %>%
    .[MT.Allele %like% "HLA", spc := "Hu"]

  dtls <- dti %>% split(by = "spc")

dtloo <- lapply(dtls, function(dti){

dtlo <- lapply(8:15, function(nmerl){

    dti <- dti[nchar(WT.Peptide) == nmerl & nchar(MT.Peptide) == nmerl]

    db <- dti[, spc %>% unique]

    dti[, spc := NULL]

    if (db == "Ms") db <- "~/antigen.garnish/Mu_iedb.fasta"
    if (db == "Hu") db <- "~/antigen.garnish/iedb.bdb"

    if (nrow(dti) == 0) {
      warning(paste("No ", nmerl, "mers compatible for fitness modeling.", sep = ""))
      return(NULL)
    }

    # drop 8mers because python scripts check anchor residues at 2 and 9 and want to avoid errors
    if (nmerl == 8) return(NULL)

    dti <- dti[, ID := link_uuid, by = 1:nrow(dti)]

    # make compatible IDs
      v <- dti[, ID %>% unique]
      names(v) <- seq_along(v)

      dti[, ID := ID %>% (function(id){
          id <- names(v[which(id == v)])
        }), by = "ID"]

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

      dn <- paste("fitness_model_", nmerl, sep = "")

      if (!dir.exists(dn)) dir.create(dn, showWarnings = FALSE)

      lapply(dtl %>% seq_along, function(i){
        dti <- data.table::copy(dtl[[i]])
        dti[, label_wt := paste(Sample, "WT", ID, MUTATION_ID, sep = "|")]
        dti[, label_mut := paste(Sample, "MUT", ID, MUTATION_ID, sep = "|")]
        aa <- dti[, WT.Peptide]
        names(aa) <- dti[, label_wt]
        aa_2 <- dti[, MT.Peptide]
        names(aa_2) <- dti[, label_mut]
        aa <- Biostrings::AAStringSet(c(aa, aa_2), use.names = TRUE)
        filename <- paste0(dn, "/neoantigens_", dti[, Sample %>% unique], ".fasta")
        aa %>% Biostrings::writeXStringSet(filename)

    # run blastp on fasta

      system(
        paste(
      "blastp -query", filename, "-db", db, "-num_threads", parallel::detectCores(), "-outfmt 5 -evalue 100000000 -gapopen 11 -gapextend 1 >",
      filename %>% stringr::str_replace("\\.fasta", replacement = "_iedb.xml"),
        sep = " ")
      )

    })

  # run pipeline

    py_path <- system.file(package = "antigen.garnish") %>% file.path(., "extdata/src/main.py")

    on <- paste(nmerl, "_neoantigens_fitness_model_output.txt", sep = "")

    system(
    paste("python",
          py_path,
          "neoepitopes.txt",
          dn,
          a,
          k,
          on,
          nmerl
    )
    )

  # curate output
    if (!file.exists(on)) {
    warning(paste("garnish_fitness did not return output for nmers of length", nmerl, sep = " "))
    return(NULL)
    }

    dto <- on %>% data.table::fread

    dto %>% data.table::setnames(c("Mutation", "Sample"),
                               c("var_uuid", "sample_id"))

    dto[, R := max(R), by = c("MutantPeptide")]

    dto %<>% melt(measure.vars = c("WildtypePeptide", "MutantPeptide"),
                value.name = "nmer")

    dt <- merge(dt[nchar(nmer) == nmerl], dto[Excluded == FALSE, .SD %>% unique,
          .SDcols = c("var_uuid",
                      "sample_id",
                      "ResidueChangeClass",
                      "R",
                      "NeoantigenRecognitionPotential",
                      "nmer")],
                      all.x = TRUE,
                      by = c("var_uuid",
                             "sample_id",
                             "nmer"))

    return(dt)

  }) %>% data.table::rbindlist(fill = TRUE) %>%
          list(., dt[nchar(nmer) == 8]) %>%
            data.table::rbindlist(fill = TRUE)

  return(dtlo)

  })

  if (length(dtls) == 2) dtloo %<>% data.table::rbindlist

  if (length(dtls) == 1) dtloo <- dtloo[[1]]

  return(dtloo)

}
