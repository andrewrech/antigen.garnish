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
#' * **NeoantigenRecognitionPotential**: Neoantigen Recognition Potential calculated by the model.
#'
#' @export garnish_fitness
#' @md

garnish_fitness <- function(dt,
                            a = 26,
                            k = 4.86936){

  # magrittr version check, this will not hide the error, only the NULL return on successful exit
  invisible(check_dep_versions())

    on.exit({
          message("Removing temporary files")
          try(
          list.files(pattern = "neoepitopes.txt") %>% file.remove, silent = TRUE)
          try(
          unlink(list.files(pattern = "Lukza_model_[0-9]+"), recursive = TRUE, force = TRUE), silent = TRUE)
                              })

  if (identical(Sys.getenv("TESTTHAT"), "true")) setwd("~")

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

	# hold out non-9mers
	dt_holdout <- dt[nchar(nmer) != 9]

	# perform Lukza modeling

    dti <- dti[nchar(WT.Peptide) == 9 & nchar(MT.Peptide) == 9]

    if (nrow(dti) == 0) {
      warning(paste("No ", 9, "mers compatible for Lukza et al. Nature 2017 fitness modeling.", sep = ""))
      return(dt_holdout)
    }

    db <- dti[, spc %>% unique]

    dti[, spc := NULL]

    if (db == "Ms") db <- "antigen.garnish/Mu_iedb.fasta"
    if (db == "Hu") db <- "antigen.garnish/iedb.bdb"

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

      dn <- paste("Lukza_model_", 9, sep = "")

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
    # https://www.ncbi.nlm.nih.gov/books/NBK279684/
    # flags here taken from Lukza et al.:
    # -task blastp-short optimized blast for <30 AA, uses larger word sizes
    # -matrix use BLOSUM62 sub substitution matrix
    # -evalue expect value for saving hits
    # -gapopen, -gapextend, numeric cost to a gapped alignment and
    # -outfmt, out put a csv with colums, seqids for query and database seuqnence, start and end of sequence match,
    # length of overlap, number of mismatches, percent identical, expected value, bitscore

      system(
        paste(
      "blastp -query", filename, "-db", db, "-num_threads", parallel::detectCores(), "-outfmt 5 -evalue 100000000 -gapopen 11 -gapextend 1 >",
      filename %>% stringr::str_replace("\\.fasta", replacement = "_iedb.xml"),
        sep = " ")
      )

    })

  # run pipeline

    py_path <- system.file(package = "antigen.garnish") %>% file.path(., "extdata/src/main.py")

    on <- paste(9, "_neoantigens_Lukza_model_output.txt", sep = "")

    system(
    paste("python",
          py_path,
          "neoepitopes.txt",
          dn,
          a,
          k,
          on,
         9
    )
    )

  # curate output
    if (!file.exists(on)) {
    warning(paste("garnish_fitness did not return output for nmers of length", 9, sep = " "))
    return(dt_holdout)
    }

    dto <- on %>% data.table::fread

    dto %>% data.table::setnames(c("Mutation", "Sample"),
                               c("var_uuid", "sample_id"))

    dto[, R := max(R), by = c("MutantPeptide")]

    dto %<>% melt(measure.vars = c("WildtypePeptide", "MutantPeptide"),
                value.name = "nmer")

    dt <- merge(dt[nchar(nmer) == 9], dto[Excluded == FALSE, .SD %>% unique,
          .SDcols = c("var_uuid",
                      "sample_id",
                      "NeoantigenRecognitionPotential",
                      "nmer")],
                      all.x = TRUE,
                      by = c("var_uuid",
                             "sample_id",
                             "nmer"))

  # re-combine 9-mers and non-9-mers
  dtlo <- list(dt, dt_holdout) %>%
       data.table::rbindlist(fill = TRUE)

  return(dtlo)

  })

  if (length(dtls) == 2) dtloo %<>% data.table::rbindlist %>% unique
  if (length(dtls) == 1) dtloo <- dtloo[[1]] %>% unique

  return(dtloo)

}


## ---- clone_wars
#' Internal function to integrate clonality data into final garnish_score.
#'
#' Integrates clonality input to create a summary metric of tumor fitness similar to [Luksza et al. *Nature* 2017](https://www.ncbi.nlm.nih.gov/pubmed/29132144).
#'
#' @import mclust
#'
#' @return A data table with added fitness model parameter columns:
#' * **clone_id**: The rank of the clone containing the variant in that sample, with the first being the largest fraction of the tumor.
#' * **clone_prop**: The estimated clustered mean value for the proportion of the tumor composed of that clone.
#' * **garnish_score**: the summary parameter of immunogenicity at the sample level, summed across dominant neoepitopes of each clone.
#'
#' @export clone_wars
#' @md

clone_wars <- function(dt){

  if (!"CELLFRACTION" %chin% names(dt) & !"AF" %chin% names(dt)){

    message("No CELLFRACTION or AF column found, returning dt without computing garnish_score from clonality.")
    return(dt)

  }

  if ("AF" %chin% names(dt)) col <- "AF"
  if ("CELLFRACTION" %chin% names(dt)) col <- "CELLFRACTION"

    b <- data.table::copy(dt)

    b %>% setnames(col, "cf")

    match_clone <- function(cf, v){

      dt <- lapply(cf, function(x){

        abs <- abs(x - v)

        a <- v[which(abs == min(abs))]
        b <- names(v)[which(abs == min(abs))]

        return(data.table(clone_prop = a, clone_id = b))

      }) %>% data.table::rbindlist

      return(dt)

    }

    cdt <- lapply(b[, sample_id %>% unique], function(s){

      dt <- b[!is.na(cf) & sample_id == s, cf, by = "var_uuid"] %>% unique

      if (nrow(dt) == 0) return(NULL)
      if (nrow(dt)== 1) return(dt[,  clone_id := 1] %>% .[, clone_prop := cf])

      x <-  mclust::Mclust(dt[, cf], verbose = FALSE)

      vect <- x$parameters$mean

      vect %<>% sort(decreasing = TRUE)

      names(vect) <- 1:length(vect) %>% as.character

      dt[, c("clone_prop", "clone_id") := match_clone(cf, v = vect)]

    }) %>% data.table::rbindlist(use.names = TRUE)

    cdt %>% data.table::setnames("cf", col)

    dt <- merge(dt, cdt, all.x = TRUE, by = c("var_uuid", col))

    # if using AF as surrogate clonality, recalculate allele fractions into cell population proportions
    if (col == "AF"){

      message("Garnish_score uses clonality to generate a sample-level metric for immune fitness.
      The provided allele fractions are being substituted, this is, however, not validated.  See ?garnish_predictions and use with care.")

      a <- dt[!is.na(clone_prop), clone_prop %>% unique, by = c("sample_id", "var_uuid")]

      # assume even ploidy and monophyletic tree so max AF is assumed to be 100% cell fraction and others are subclonal (assumes no pressure to lose mutation and max AF = tumor sample purity).
      # teleologically, this assumption is really only good for a relatively genetically stable clonal cell line taken off a dish and sequenced.
      # can't see this being at all a large task, can set max cores and not worry about memory.

      a[, clone_prop := V1 / max(V1, na.rm = TRUE), by = "sample_id"]

      # clear column derived from allele fractions in dt
      dt[, clone_prop := NULL]

      dt <- merge(dt, a[, .SD %>% unique, .SDcol = c("sample_id", "var_uuid", "clone_prop")],
      all.x = TRUE, by = c("sample_id", "var_uuid"))

    }

  dt[!is.na(clone_prop), efit := exp(fitness_score %>% max(na.rm = TRUE)) * unique(clone_prop), by = c("sample_id", "clone_id")]

  dt[!is.na(efit), garnish_score := efit %>% unique %>% sum(na.rm = TRUE), by = "sample_id"]

  dt[, efit := NULL]

  return(dt)

}
