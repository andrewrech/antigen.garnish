## ---- garnish_fitness
#' Generate neoepitope fitness model
#'
#' Implements the neoepitope fitness model of [Luksza et al. *Nature* 2017](https://www.ncbi.nlm.nih.gov/pubmed/29132144).
#'
#' @param dt Data table output from `garnish_affinity`.
#' @param a Numeric fitness model parameter. Binding curve horizontal displacement used to determine TCR recognition probability of a peptide compared to an IEDB near match.
#' @param k Numeric fitness model parameter. Steepness of the binding curve at `a`.
#'
#' @return A data table with added fitness model parameter columns:
#' * **NeoantigenRecognitionPotential**: Neoantigen recognition potential calculated by the model.
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


## ---- garnish_clonality
#' Internal function to integrate clonality data into final `garnish_score` parameter.
#'
#' Integrates clonality input for creation of a summary metric of tumor fitness similar to [Luksza et al. *Nature* 2017](https://www.ncbi.nlm.nih.gov/pubmed/29132144).
#'
#' @param dt Data table. Passed internally from garnish_affinity, requires allelic_fraction or cell_fraction columns.
#'
#' @return A data table with added fitness model parameter columns:
#' * **clone_id**: rank of the clone containing the variant (highest equals larger tumor fraction).
#' * **cl_proportion**: estimated mean tumor fraction containing the clone.
#' * **garnish_score**: the summary parameter of immunogenicity at the sample level, summed across top neoepitopes of each clone.
#'
#' @export garnish_clonality
#' @md

garnish_clonality <- function(dt){

  if (!"cell_fraction" %chin% names(dt) & !"allelic_fraction" %chin% names(dt)){

    warnings("No cell_fraction or allelic_fraction column found. Returning dt without computing garnish_score from clonality.")
    return(dt)
  }

  if ("allelic_fraction" %chin% names(dt))
  	col <- "allelic_fraction"

  # prefer cell_fraction if available
  if ("cell_fraction" %chin% names(dt))
  	col <- "cell_fraction"


    b <- data.table::copy(dt)

    b %>% setnames(col, "cf")

    match_clone <- function(cf, v){

      dt <- lapply(cf, function(x){

        abs <- abs(x - v)

        a <- v[which(abs == min(abs))]
        b <- names(v)[which(abs == min(abs))]

        return(data.table(cl_proportion = a, clone_id = b))

      }) %>% data.table::rbindlist

      return(dt)

    }

    cdt <- lapply(b[, sample_id %>% unique], function(s){

      dt <- b[!is.na(cf) & sample_id == s, cf, by = "var_uuid"] %>% unique

      if (nrow(dt) == 0)
      	return(NULL)

      if (nrow(dt)== 1)
      	return(dt[,  clone_id := 1] %>% .[, cl_proportion := cf])

      x <-  mclust::Mclust(dt[, cf], verbose = FALSE)

      vect <- x$parameters$mean

      vect %<>% sort(decreasing = TRUE)

      names(vect) <- 1:length(vect) %>% as.character

      dt[, c("cl_proportion", "clone_id") := match_clone(cf, v = vect)]

    }) %>% data.table::rbindlist(use.names = TRUE)

    cdt %>% data.table::setnames("cf", col)

    dt <- merge(dt, cdt, all.x = TRUE, by = c("var_uuid", col))

    # remove cl_proportion from wt peptide rows because its meaningless and refers to matched mutant nmer
    dt[pep_type == "wt", cl_proportion := as.numeric(NA)]

    # if using allelic_fraction as surrogate clonality, recalculate allele fractions into cell population proportions

    if (col == "allelic_fraction"){

      a <- dt[!is.na(cl_proportion) & !is.na(allelic_fraction) & pep_type != "wt", cl_proportion %>% unique, by = c("sample_id", "var_uuid", "pep_type")]

      # determine maximum allelic fraction from most prevalent SNVs, use top decile from ecdf

      ecdf_wrap <- function(v){

        return(v %>% stats::ecdf %>% stats::quantile(0.9))

      }

      a[, ecdf := ecdf_wrap(V1), by = "sample_id"]

      a[, cl_proportion := V1 / ecdf, by = "sample_id"]

      # clear column derived from allele fractions in dt
      dt[, cl_proportion := NULL]

      dt <- merge(dt, a[, .SD %>% unique, .SDcol = c("sample_id", "var_uuid", "cl_proportion", "pep_type")],
      all.x = TRUE, by = c("sample_id", "var_uuid", "pep_type"))

    }

  dt[!is.na(cl_proportion), efit := exp(fitness_score %>% max(na.rm = TRUE)) * cl_proportion, by = c("sample_id", "clone_id")]

  dt[!is.na(efit), garnish_score := efit %>% unique %>% sum(na.rm = TRUE), by = "sample_id"]

  dt[, efit := NULL]

  return(dt)

}
