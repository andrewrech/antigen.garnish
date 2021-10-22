
#' Return foreignness_scores for a vector of peptides.
#'
#' @param v Character. Vector of nmers.
#' @param db Character. One of c("mouse", "human")
#'
#' @return Data table of nmers and corresponding foreignness_score values.
#'
#' @export foreignness_score
#' @md

foreignness_score <- function(v, db) {
  if (!db %chin% c("mouse", "human")) stop("db must be \"human\" or \"mouse\"")

  on.exit({
    message("Removing temporary fasta files.")
    try(
      list.files(pattern = "foreignness_score|blastp") %>% file.remove()
    )
  })

  if (suppressWarnings(system("which blastp 2> /dev/null", intern = TRUE)) %>%
    length() == 0) {
    warning("Skipping BLAST because ncbiblast+ is not in PATH")
    return(data.table::data.table(nmer = v))
  }

  message("Generating FASTA to query.")

  names(v) <- 1:length(v) %>% as.character()

  sdt <- v %>%
    data.table::as.data.table() %>%
    .[, nmer_id := names(v)]

  AA <- Biostrings::AAStringSet(v, use.names = TRUE)
  Biostrings::writeXStringSet(AA, file = "foreignness_score_fasta.fa", format = "fasta")

  # run blastp-short for iedb matches
  # https://www.ncbi.nlm.nih.gov/books/NBK279684/
  # flags here taken from Lukza et al.:
  # -task blastp-short optimized blast for <30 AA, uses larger word sizes
  # -matrix use BLOSUM62 sub substitution matrix
  # -evalue expect value for saving hits
  # -gapopen, -gapextend, numeric cost to a gapped alignment and
  # -outfmt, out put a csv with colums, seqids for query and database seuqnence, start and end of sequence match,
  # length of overlap, number of mismatches, percent identical, expected value, bitscore

  message("Running blastp for homology to IEDB antigens.")

  ag_dir <- Sys.getenv("AG_DATA_DIR")

  if (db == "mouse") db <- paste0(ag_dir, "/Mu_iedb.fasta.pin")

  if (db == "human") db <- paste0(ag_dir, "/iedb.bdb.pin")

  if (!file.exists(db)) {
    warn <- paste(
      "Skipping foreignness_score because BLAST database cannot be found.",
      "To set a custom path to the antigen.garnish data folder",
      "set environomental variable AG_DATA_DIR from the shell",
      "or from R using Sys.setenv",
      "",
      "Re-download installation data:",
      '$ curl -fsSL "http://get.rech.io/antigen.garnish-2.3.0.tar.gz" | tar -xvz',
      "",
      "Documentation:",
      "https://neoantigens.rech.io",
      sep = "\n"
    )
    warnings(warn)
    return(data.table::data.table(nmer = v))
  }

  db <- paste("-db", db %>% stringr::str_replace("\\.pin", ""), sep = " ")

  threads <- Sys.getenv("AG_THREADS") %>% as.numeric(.)
  if (threads == "" || threads %>% is.na) {
    threads <- parallel::detectCores()
  }
  if ((threads %>% as.numeric()) %>% is.na()) {
    stop("Error reading environment variable AG_THREADS as numeric.")
  }

  system(paste0(
    "blastp -query foreignness_score_fasta.fa ", db, " -evalue 100000000 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -out blastp_iedbout.csv -num_threads ", threads,
    " -outfmt '10 qseqid sseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'"
  ))

  blastdt <- list.files(pattern = "iedbout\\.csv")

  if (length(blastdt) == 0) {
    message("No blast output against IEDB returned.")
    return(data.table::data.table(nmer = v))
  }

  if (all(file.info(blastdt)$size == 0)) {
    message("No IEDB matches found by blast.")
    return(data.table::data.table(nmer = v, foreignness_score = 0))
  }

  blastdt <- blastdt %>% data.table::fread()

  blastdt %>% data.table::setnames(
    names(.),
    c(
      "nmer_id",
      "IEDB_anno",
      "nmer",
      "q_start",
      "q_stop",
      "WT.peptide",
      "s_start",
      "s_end",
      "overlap_length",
      "mismatch_length",
      "pident",
      "evalue",
      "bitscore"
    )
  )

  blastdt <- blastdt[, nmer := nmer %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[, WT.peptide := WT.peptide %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[!is.na(nmer) & !is.na(WT.peptide)]

  blastdt <- blastdt[nmer %like% "^[ARNDCQEGHILKMFPSTWYV]+$" & WT.peptide %like% "^[ARNDCQEGHILKMFPSTWYV]+$"]

  if (nrow(blastdt) == 0) {
    message(paste("No IEDB matches found with cannonical AAs, can't compute foreignness score...."))
    return(data.table::data.table(nmer = v))
  }

  message("Summing IEDB local alignments...")

  blastdt[, SW := make_sw_alignment(nmer, WT.peptide)]

  modeleR <- function(als, a = 26, k = 4.86936) {
    be <- -k * (a - als)

    sumexp <- sum(exp(be))

    Zk <- 1 + sumexp
    R <- sumexp / Zk

    return(R)
  }

  blastdt[, foreignness_score := SW %>% modeleR(), by = "nmer_id"]

  # get full IEDB ref here
  fa <- Biostrings::readAAStringSet(db %>%
    stringr::str_replace("bdb$", "fasta") %>%
    stringr::str_replace("^-db\\ ", ""))
  f <- fa %>% as.character()
  names(f) <- names(fa)

  blastdt[, IEDB_anno := lapply(IEDB_anno, function(i) {
    mv <- f[which(stringr::str_detect(pattern = stringr::fixed(i), names(f)))]
    mv <- mv[which(stringr::str_detect(pattern = stringr::fixed(WT.peptide), mv))]
    return(paste(names(mv), WT.peptide, collapse = "|"))
  }), by = 1:nrow(blastdt)]

  anndt <- blastdt[, .SD %>% unique(), .SDcols = c(
    "nmer_id",
    "nmer",
    "IEDB_anno",
    "foreignness_score",
    "SW"
  )]

  blastdt <- blastdt[, .SD %>% unique(), .SDcols = c("nmer_id", "foreignness_score")]

  sdt[, nmer_id := as.character(nmer_id)]
  blastdt[, nmer_id := as.character(nmer_id)]

  sdt <- merge(sdt, blastdt, by = "nmer_id")

  sdt %>% data.table::setnames(".", "nmer")

  anndt <- anndt[, msw := max(SW), by = "nmer_id"] %>%
    .[SW == msw] %>%
    .[, .SD %>% unique(), .SDcols = c("nmer_id", "IEDB_anno")]

  # merge equally good IEDB_annos into one
  anndt[, IEDB_anno := paste(IEDB_anno %>% unique(), collapse = "|"), by = "nmer_id"]

  anndt %<>% unique

  anndt[, nmer_id := as.character(nmer_id)]

  sdt <- merge(sdt, anndt, by = "nmer_id", all.x = TRUE)

  sdt <- sdt[, .SD %>% unique(), .SDcols = c("nmer", "foreignness_score", "IEDB_anno")]

  return(sdt)
}

#' Return dissimilarity (to reference proteome) values for a vector of peptides.
#'
#' @param v Character. Vector of nmers.
#' @param db Character. One of c("mouse", "human").
#' @param kval Numeric. Steepness of sigmoidal curve at k. Default 4.86936, the value used in the analysis of Van Allen, Snyder, Rizvi, Riaz, and Hellmann datasets.
#' @param aval Numeric. Optionally can be "mean" to use mean alignment for nmers passed. Horizontal displacement of partition function. Default is 32, based on max_SW of 75 million 8-15mers from the five clinical datasets against human, if using max_SW, use 52. This value may not be meaningful for murine alignment so use with care.
#'
#' @return Data table of nmers and corresponding dissimilarity values (to the non-mutated proteome).
#'
#' @export dissimilarity_score
#' @md

dissimilarity_score <- function(v, db, kval = 4.86936, aval = 32) {
  if (!db %chin% c("mouse", "human")) stop("db must be \"human\" or \"mouse\"")

  on.exit({
    message("Removing temporary fasta files.")
    try(
      list.files(pattern = "dissimilarity|blastp|blastdt_[0-9]+\\.txt$") %>% file.remove()
    )
  })

  if (identical(Sys.getenv("TESTTHAT"), "true")) setwd(Sys.getenv("HOME"))

  if (suppressWarnings(system("which blastp 2> /dev/null", intern = TRUE)) %>%
    length() == 0) {
    warning("Skipping BLAST because ncbiblast+ is not in PATH")
    return(data.table::data.table(nmer = v))
  }

  # generate fastas to query

  names(v) <- 1:length(v) %>% as.character()

  sdt <- v %>%
    data.table::as.data.table() %>%
    .[, nmer_id := names(v)]

  AA <- Biostrings::AAStringSet(v, use.names = TRUE)
  Biostrings::writeXStringSet(AA, file = "dissimilarity_fasta.fa", format = "fasta")

  # run blastp for iedb matches
  # https://www.ncbi.nlm.nih.gov/books/NBK279684/
  # flags here taken from Lukza et al.:
  # -matrix use BLOSUM62 substitution matrix
  # -evalue expect value for saving hits
  # -gapopen, -gapextend, numeric cost to a gapped alignment and
  # -outfmt, out put a csv with colums, seqids for query and database seuqnence, start and end of sequence match,
  # length of overlap, number of mismatches, percent identical, expected value, bitscore
  message("Running blastp for homology to self antigens.")

  ag_dir <- Sys.getenv("AG_DATA_DIR")

  if (db == "mouse") db <- paste0(ag_dir, "/mouse.bdb.pin")

  if (db == "human") db <- paste0(ag_dir, "/human.bdb.pin")

  # check one level up for blastdb, useful for test env temp folders
  if (!file.exists(db)) db <- file.path("..", db)

  if (!file.exists(db)) {
    warn <- paste(
      "Skipping dissimilarity because BLAST database cannot be found.",
      "To set a custom path to the antigen.garnish data folder",
      "set environomental variable AG_DATA_DIR from the shell",
      "or from R using Sys.setenv",
      "",
      "Re-download installation data:",
      '$ curl -fsSL "http://get.rech.io/antigen.garnish-2.3.0.tar.gz" | tar -xvz',
      "",
      "Documentation:",
      "https://neoantigens.rech.io",
      sep = "\n"
    )
    warnings(warn)
    return(data.table::data.table(nmer = v))
  }

  if (db == "mouse") {
    message(paste(
      "Dissimilarity was modeled on 75 million missense derived",
      "peptides against the human proteome,",
      "applicability to murine data is unknown."
    ))
  }


  db <- paste("-db", db %>% stringr::str_replace("\\.pin", ""), sep = " ")

  threads <- Sys.getenv("AG_THREADS") %>% as.numeric(.)
  if (threads == "" || threads %>% is.na) {
    threads <- parallel::detectCores()
  }
  if ((threads %>% as.numeric()) %>% is.na()) {
    stop("Error reading environment variable AG_THREADS as numeric.")
  }

  system(paste0(
    "blastp -query dissimilarity_fasta.fa ", db, " -evalue 100000000 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -out blastp_self.csv -num_threads ", threads,
    " -outfmt '10 qseqid sseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'"
  ))


  blastdt <- list.files(pattern = "blastp_self\\.csv")

  if (length(blastdt) == 0) {
    warning("No blast output for dissimilarity returned.")
    return(data.table::data.table(nmer = v))
  }

  if (all(file.info(blastdt)$size == 0)) {
    message("No self-proteome matches found by blast.")
    return(data.table::data.table(nmer = v, dissimilarity = 0))
  }

  blastdt <- blastdt %>% data.table::fread()

  blastdt %>% data.table::setnames(
    names(.),
    c(
      "nmer_id",
      "self_anno",
      "nmer",
      "q_start",
      "q_stop",
      "WT.peptide",
      "s_start",
      "s_end",
      "overlap_length",
      "mismatch_length",
      "pident",
      "evalue",
      "bitscore"
    )
  )

  blastdt <- blastdt[, nmer := nmer %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[, WT.peptide := WT.peptide %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[!is.na(nmer) & !is.na(WT.peptide)]

  blastdt <- blastdt[nmer %like% "^[ARNDCQEGHILKMFPSTWYV]+$" & WT.peptide %like% "^[ARNDCQEGHILKMFPSTWYV]+$"]

  if (nrow(blastdt) == 0) {
    message(paste("No self matches found with cannonical AAs, can't compute dissimilarity...."))
    return(data.table::data.table(nmer = v))
  }

  message("Calculating dissimilarity...")

  # this is memory intense, lets split it up and stitch it back
  suppressWarnings(
    blastdt <- blastdt %>% split(1:(nrow(blastdt) / 100))
  )

  blastdt <- lapply(blastdt %>% seq_along(), function(i) {
    fn <- paste("blastdt_", i, ".txt", sep = "")

    blastdt[[i]] %>% data.table::fwrite(fn, sep = "\t")

    return(fn)
  }) %>% unlist()

  lapply(blastdt %>% seq_along(), function(i) {

    # print(paste("Alignment subset", i, "of", length(blastdt)))

    b <- blastdt[i] %>% data.table::fread()

    b[, SW := make_sw_alignment(nmer, WT.peptide)]

    b %>% data.table::fwrite(blastdt[i], sep = "\t")

    return(NULL)
  })

  blastdt <- lapply(blastdt, function(f) {
    dt <- data.table::fread(f)

    file.remove(f)

    return(dt)
  }) %>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

  message("Running partition function...")

  modeleR <- function(als, a = aval, k = kval) {
    be <- -k * (a - als)

    sumexp <- sum(exp(be))

    Zk <- 1 + sumexp
    R <- sumexp / Zk

    R <- 1 - R

    return(R)
  }

  blastdt[, dissimilarity := SW %>% modeleR(), by = "nmer_id"]

  blastdt <- blastdt[, .SD %>% unique(), .SDcols = c("nmer_id", "dissimilarity")]

  sdt[, nmer_id := as.character(nmer_id)]
  blastdt[, nmer_id := as.character(nmer_id)]

  sdt <- merge(sdt, blastdt, by = "nmer_id")

  sdt %>% data.table::setnames(".", "nmer")

  return(sdt[, .SD %>% unique(), .SDcols = c("nmer", "dissimilarity")])
}


#' Internal function to categorize mutant and wild-type peptides by similarity using `BLAST` to calculate neoantigen amplitude and homology to IEDB antigens.
#'
#' @param dti Data table of nmers.
#'
#' @noRd

make_BLAST_uuid <- function(dti) {
  on.exit({
    message("Removing temporary fasta files.")
    try(
      list.files(pattern = "(Ms|Hu)_nmer_fasta") %>% file.remove()
    )
    try(
      list.files(pattern = "dai_blastdt_[0-9]+") %>% file.remove()
    )
    try(blast_out %>% file.remove())
  })

  if (suppressWarnings(system("which blastp 2> /dev/null", intern = TRUE)) %>%
    length() == 0) {
    warning("Skipping BLAST because ncbiblast+ is not in PATH")
    return(dti)
  }

  # check for data directory
  check_pred_tools()

  dt <- dti[pep_type != "wt" & !is.na(pep_type)]

  # blast first to get pairs for non-mutnfs peptides then run nature paper package

  dt[MHC %like% "H-2" | MHC == "all_mouse", spc := "Ms"] %>%
    .[MHC %like% "HLA" | MHC == "all_human", spc := "Hu"]

  # generate fastas to query

  lapply(dt[, spc %>% unique()], function(s) {
    dt <- dt[spc == s]

    fa_v <- dt[, .SD %>% unique(), .SDcols = c("nmer", "nmer_uuid")] %>%
      .[order(nmer)] %>%
      .[, nmer]

    names(fa_v) <- dt[, .SD %>%
      unique(), .SDcols = c("nmer", "nmer_uuid")] %>%
      .[order(nmer)] %>%
      .[, nmer_uuid]

    AA <- Biostrings::AAStringSet(fa_v, use.names = TRUE)
    Biostrings::writeXStringSet(AA, file = paste(s, "_ag_nmer_fasta.fa", sep = ""), format = "fasta")
  })

  # run blastp-short for near matches
  # https://www.ncbi.nlm.nih.gov/books/NBK279684/
  # flags here indicate:
  # -task blastp-short optimized blast for <30 AA, uses larger word sizes
  # -outfmt, out put a csv with colums, seqids for query and database seuqnence, start and end of sequence match,
  # length of overlap, number of mismatches, percent identical, expected value, bitscore
  # PAM30 is default substitution matrix here

  message("Running blastp-short to find close matches for differential agretopicity calculation.")

  blast_out <- paste(dt[, spc %>% unique()],
    uuid::UUIDgenerate() %>% substr(1, 18) %>% stringr::str_replace_all("-", ""),
    "_blastpout.csv",
    sep = ""
  )

  threads <- Sys.getenv("AG_THREADS") %>% as.numeric(.)
  if (threads == "" || threads %>% is.na) {
    threads <- parallel::detectCores()
  }
  if ((threads %>% as.numeric()) %>% is.na()) {
    stop("Error reading environment variable AG_THREADS as numeric.")
  }


  if (file.exists("Ms_ag_nmer_fasta.fa")) {
    cmd <- paste0(
      "blastp -query Ms_ag_nmer_fasta.fa -task blastp-short -db ",
      Sys.getenv("AG_DATA_DIR"),
      "/mouse.bdb -out ",
      blast_out[blast_out %like% "Ms"],
      " -num_threads ", threads,
      " -outfmt '10 qseqid sseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'"
    )

    cmd %>% message()

    result <- cmd %>% system()

    if (!result == 0) {
      warnings("Blastp-short to find close matches for differential agretopicity calculation had a non-zero exit status")
    }
  }

  if (file.exists("Hu_ag_nmer_fasta.fa")) {
    cmd <- paste0(
      "blastp -query Hu_ag_nmer_fasta.fa -task blastp-short -db ",
      Sys.getenv("AG_DATA_DIR"),
      "/human.bdb -out ",
      blast_out[blast_out %like% "Hu"],
      " -num_threads ", threads,
      " -outfmt '10 qseqid sseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'"
    )

    cmd %>% message()

    result <- cmd %>% system()

    if (!result == 0) {
      warnings("Blastp-short to find close matches for differential agretopicity calculation had a non-zero exit status")
    }
  }


  # read in CSV

  blastdt <- blast_out

  if (length(blastdt) == 0 ||
    all(file.info(blastdt)$size == 0)) {
    message("No WT similarity matches found by blast.")
    return(dti)
  }

  if (any(file.info(blastdt)$size == 0)) {
    blastdt <- blastdt %>%
      file.info() %>%
      data.table::as.data.table(keep.rownames = TRUE) %>%
      .[size != 0, rn]
  }


  blastdt <- lapply(blastdt, data.table::fread) %>%
    data.table::rbindlist(use.names = TRUE) %>%
    data.table::setnames(
      names(.),
      c(
        "nmer_uuid",
        "ensembl_prot_id",
        "nmer",
        "q_start",
        "q_stop",
        "WT.peptide",
        "s_start",
        "s_end",
        "overlap_length",
        "mismatch_length",
        "pident",
        "evalue",
        "bitscore"
      )
    )

  if (nrow(blastdt) == 0) {
    message("No WT similarity matches found by blast.")
    return(dti)
  }

  blastdt <- blastdt[, nmer := nmer %>%
    stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[, WT.peptide := WT.peptide %>%
      stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[nchar(nmer) != 0 & nchar(WT.peptide) != 0] %>%
    .[nchar(nmer) == nchar(WT.peptide) & mismatch_length == 1] %>%
    .[!is.na(nmer) & !is.na(WT.peptide)]

  # remove uncertain AA calls/whitelist cannonical AA
  blastdt <- blastdt[nmer %like% "^[ARNDCQEGHILKMFPSTWYV]+$" & WT.peptide %like% "^[ARNDCQEGHILKMFPSTWYV]+$"]

  if (nrow(blastdt) == 0) {
    message("No WT similarity matches found by blast.")
    return(dti)
  }

  message("Calculating local alignment to WT peptides for proteome-wide differential agretopicity predictions.")

  # this is memory intense, lets split it up and stitch it back
  suppressWarnings(
    blastdt <- blastdt %>% split(1:(nrow(blastdt) / 100))
  )

  blastdt <- lapply(blastdt %>% seq_along(), function(i) {
    fn <- paste("dai_blastdt_", i, ".txt", sep = "")

    blastdt[[i]] %>% data.table::fwrite(fn, sep = "\t")

    return(fn)
  }) %>% unlist()

  lapply(blastdt %>% seq_along(), function(i) {
    print(paste("Alignment subset", i, "of", length(blastdt)))

    b <- blastdt[i] %>% data.table::fread()

    b[, SW := make_sw_alignment(nmer, WT.peptide)]

    b %>% data.table::fwrite(blastdt[i], sep = "\t")

    return(NULL)
  })

  blastdt <- lapply(blastdt, function(f) {
    dt <- data.table::fread(f)

    file.remove(f)

    return(dt)
  }) %>% data.table::rbindlist(use.names = TRUE)

  blastdt[, highest := max(SW), by = "nmer_uuid"]

  # keep highest local alignment
  blastdt <- blastdt[highest == SW] %>%
    .[, highest := NULL]

  # keep longest alignment
  blastdt <- blastdt[, best := max(overlap_length), by = "nmer_uuid"] %>%
    .[best == overlap_length] %>%
    .[, best := NULL]

  # dedupe but retain multiple equally good SW scoring matches by sequence (not by match source)
  blastdt %>% data.table::setkey(nmer_uuid, WT.peptide)
  blastdt %<>% unique(by = c("nmer_uuid", "WT.peptide"))

  blastdt <- blastdt[, .SD %>% unique(), .SDcols = c("nmer_uuid", "nmer", "WT.peptide")]
  blastdt[, blast_uuid := uuid::UUIDgenerate(), by = c("nmer_uuid", "WT.peptide")]

  dti <- merge(dti,
    blastdt[, .SD, .SDcols = c("nmer_uuid", "nmer", "blast_uuid")],
    by = c("nmer_uuid", "nmer"), all.x = TRUE
  )

  # to add WT.peptide back to table need nmer, nmer_i, nmer_l (nchar(nmer)), var_uuid, pep_type

  if (!"var_uuid" %chin% names(dti)) dti[, var_uuid := NA %>% as.character()]
  if (!"effect_type" %chin% names(dti)) dti[, effect_type := NA %>% as.character()]

  vdt <- dti[, .SD %>% unique(), .SDcols = c("nmer_uuid", "nmer", "nmer_i", "nmer_l", "var_uuid", "sample_id", "effect_type", "MHC")]

  vdt <- merge(vdt, blastdt, by = c("nmer", "nmer_uuid"))

  vdt %>%
    .[, nmer := NULL] %>%
    .[, nmer_uuid := NULL]

  vdt <- vdt[, nmer_l := nchar(WT.peptide)] %>%
    data.table::setnames("WT.peptide", "nmer") %>%
    .[, nmer_uuid := uuid::UUIDgenerate(), by = c("nmer", "var_uuid")] %>%
    unique() %>%
    .[, pep_type := "wt"]

  dto <- data.table::rbindlist(list(dti, vdt), fill = TRUE, use.names = TRUE)

  return(dto)
}


#' Internal function to pair peptides from missense sites by a common UUID for DAI calculations
#'
#' @param dt Data table of nmers.
#'
#' @noRd
#' @md

make_DAI_uuid <- function(dt) {
  if (!c(
    "pep_type", "nmer", "nmer_i",
    "nmer_l", "var_uuid", "frameshift",
    "pep_mut", "pep_wt"
  ) %chin%
    (dt %>% names()) %>% any()) {
    stop("dt is missing columns")
  }

  daidt <- dt %>%
    # explicitly select missense
    .[frameshift == FALSE &
      (pep_mut %>% nchar()) == (pep_wt %>% nchar())] %>%
    # recombine table, creating pairs
    # merge vs. sort and bind so edges cases
    # result in lost peptides vs. total mis-
    # alignment
    {
      merge(
        .[pep_type == "wt", .SD,
          .SDcols = c(
            "nmer", "nmer_i",
            "nmer_l", "var_uuid"
          )
        ] %>%
          data.table::setnames("nmer", "wt_nmer"),

        .[pep_type == "mutnfs", .SD,
          .SDcols = c(
            "nmer", "nmer_i",
            "nmer_l", "var_uuid"
          )
        ] %>%
          data.table::setnames("nmer", "mtnfs_nmer"),
        by = c("var_uuid", "nmer_i", "nmer_l")
      )
    }

  daidt[, dai_uuid := uuid::UUIDgenerate(use.time = FALSE, n = .N)]

  daidt %<>%
    # bind back into one table
    {
      rbindlist(list(
        .[, .SD, .SDcols = c(
          "wt_nmer", "nmer_i",
          "nmer_l", "var_uuid", "dai_uuid"
        )] %>%
          data.table::setnames("wt_nmer", "nmer"),
        .[, .SD, .SDcols = c(
          "mtnfs_nmer", "nmer_i",
          "nmer_l", "var_uuid", "dai_uuid"
        )] %>%
          data.table::setnames("mtnfs_nmer", "nmer")
      ))
    }

  # merge back together
  dt %<>% merge(daidt,
    by = c(
      "nmer",
      "nmer_i",
      "nmer_l",
      "var_uuid"
    ),
    all.x = TRUE
  )
  return(dt)
}


#' Internal function to merge input data table and prediction results
#'
#' @param l Output list from run_netMHC
#' @param dt Input data table.
#' @noRd
#' @md

merge_predictions <- function(l, dt) {
  message("Merging output.")

  # merge netMHC by program type

  progl <- lapply(l %>% seq_along(), function(dti) {
    l[[dti]]$command[1] %>% stringr::str_extract("net[A-Za-z]+")
  })

  # merge netMHC output

  for (ptype in (progl %>% unique() %>% unlist())) {
    dt <- merge(dt, l[(progl == ptype) %>% which()] %>%
      data.table::rbindlist(),
    by = c("nmer", ptype), all.x = TRUE,
    allow.cartesian = TRUE
    )
  }

  dt %<>% unique

  message("Reading mhcflurry output.")

  f_flurry <- list.files(pattern = "mhcflurry_output.*csv")
  if (f_flurry %>% length() > 0) {
    fdt <- lapply(f_flurry, function(x) {
      dt <- suppressWarnings(data.table::fread(x))

      if (nrow(dt) == 0) {
        return(NULL)
      }

      return(dt)
    }) %>%
      data.table::rbindlist() %>%
      data.table::setnames(c("allele", "peptide"), c("MHC", "nmer"))
    dt <- merge(dt, fdt %>% unique(), by = c("nmer", "MHC"), all.x = TRUE)
    dt %<>% unique
  }

  message("Calculating netMHC consensus score.")
  for (col in (dt %>% names() %include% "aff|[Rr]ank|Ensemble_score")) {
    suppressWarnings({
      set(dt, j = col, value = dt[, get(col) %>% as.numeric()])
    })
  }

  message("Calculating overall consensus affinity score.")

  # get vector of netMHC scores
  cols <- dt %>% names() %include% c("affinity\\(nM\\)")

  # only calculate best_netMHC if 2 or more scores exist

  if (length(cols) < 2) {
    dt[, best_netMHC := get(cols)]
  }

  if (length(cols) >= 2) {
    dtm <- dt[, .SD, .SDcols = c("nmer", "MHC", cols)] %>%
      melt(id.vars = c("nmer", "MHC")) %>%
      # order affinity predictions by program preference
      .[, variable := variable %>% factor(levels = cols)] %>%
      # key table so first non-NA value is the preferred program
      data.table::setkey(nmer, MHC, variable) %>%
      .[, .(
        best_netMHC =
          # define Consensus_score
        value %>%
          na.omit() %>%
          .[1]
      ), by = c("nmer", "MHC")] %>%
      .[!best_netMHC %>% is.na()]

    dt %<>% merge(dtm, by = c("nmer", "MHC"), all = TRUE)
  }

  # merge back

  # take average of mhcflurry best available netMHC tool

  cols <- dt %>% names() %include% "(best_netMHC)|(mhcflurry_prediction$)|(mhcflurry_affinity$)"

  if (length(cols) < 2) {
    dt[, Ensemble_score := get(cols)]
  }

  if (length(cols) >= 2) {
    dt[, Ensemble_score := mean(as.numeric(.SD), na.rm = TRUE),
      by = 1:nrow(dt), .SDcols = cols
    ]
  }

  message("Calculating differential agretopicity.")

  dt[, DAI := NA %>% as.numeric()]

  data.table::setkey(dt, pep_type, dai_uuid)

  # dai_uuid is always length 2
  # calculate DAI by dividing
  # correct WT affinity per Lukza et al., "* (1 / (1 + (0.0003 * Ensemble_score[2])))".
  # only affects large WT kD that would otherwise overestimate amplitude value (affinity prediction algorithms are right skewed.)

  dt[!dai_uuid %>% is.na(),
    DAI := (Ensemble_score[2] /
      Ensemble_score[1]) *
      (1 / (1 + (0.0003 * Ensemble_score[2]))),
    by = c("dai_uuid", "MHC")
  ]

  message("Determining peptide pairs for BLAST.")

  if ("blast_uuid" %chin% names(dt)) {

    # keep blast match that will give most conservative BLAST_A value
    # suppress empty vector warning returning Inf
    suppressWarnings({
      dt[!is.na(blast_uuid) & pep_type == "wt",
        match := Ensemble_score %>% as.numeric() %>% min(na.rm = TRUE),
        by = c("blast_uuid", "MHC")
      ] %>%
        .[Ensemble_score == match, match := 0]
    })

    dt <- dt[is.na(match) | match == 0]

    dt[, match := NULL]

    data.table::setkey(dt, pep_type, blast_uuid)

    dt[!blast_uuid %>% is.na(),
      BLAST_A := (Ensemble_score[2] /
        Ensemble_score[1]) *
        (1 / (1 + (0.0003 * Ensemble_score[2]))),
      by = c("blast_uuid", "MHC")
    ]
  }

  if ("forn_uuid" %chin% names(dt)) {

    # keep blast match that will give most conservative IEDB_A value
    # suppress empty vector warning returning Inf
    suppressWarnings({
      dt[!is.na(forn_uuid) & effect_type == "IEDB_source",
        match := Ensemble_score %>% as.numeric() %>% min(na.rm = TRUE),
        by = c("forn_uuid", "MHC")
      ] %>%
        .[Ensemble_score == match, match := 0]
    })

    dt <- dt[is.na(match) | match == 0]

    dt[, match := NULL]

    data.table::setkey(dt, pep_type, forn_uuid)

    dt[!forn_uuid %>% is.na(),
      IEDB_A := (Ensemble_score[2] /
        Ensemble_score[1]) *
        (1 / (1 + (0.0003 * Ensemble_score[2]))),
      by = c("forn_uuid", "MHC")
    ]
  }

  return(dt)
}


#' Internal function to create commands for neoantigen prediction.
#'
#' @param dt Data.table of predictions to run.
#' @noRd

get_pred_commands <- function(dt) {
  if (!c("nmer", "MHC", "nmer_l") %chin%
    (dt %>% names()) %>% any()) {
    stop("dt is missing columns")
  }

  # replace all with all MHC types
  dt[MHC == "all_human", MHC :=
    system.file("extdata",
      "all_alleles.txt",
      package = "antigen.garnish"
    ) %>%
    data.table::fread(header = FALSE, sep = "\t") %>%
    .[V1 %like% "HLA"] %>%
    .$V1 %>%
    paste(collapse = " ")]

  dt[MHC == "all_mouse", MHC :=
    system.file("extdata",
      "all_alleles.txt",
      package = "antigen.garnish"
    ) %>%
    data.table::fread(header = FALSE, sep = "\t") %>%
    .[V1 %like% "H-2"] %>%
    .$V1 %>%
    paste(collapse = " ")]

  if (dt[, MHC %>% unique()] %>%
    stringr::str_detect(" ") %>% any()) {
    dt %<>%
      tidyr::separate_rows("MHC", sep = " ")
  }
  dt <- data.table::copy(dt) %>% data.table::as.data.table(.)

  dt[, class := "none"]
  dt[MHC %>% stringr::str_detect("(HLA-[ABC]\\*)|(H-2-[A-Z][a-z])"), class := "I"]
  dt[class != "I", class := "II"]

  # get available MHC alleles for predictions

  alleles <- data.table::rbindlist(
    list(
      system.file("extdata",
        "netMHC_alleles.txt",
        package = "antigen.garnish"
      ) %>%
        data.table::fread(header = FALSE, sep = "\t") %>%
        data.table::setnames("V1", "allele") %>%
        .[, type := "netMHC"],
      system.file("extdata",
        "netMHCpan_alleles.txt",
        package = "antigen.garnish"
      ) %>%
        data.table::fread(header = FALSE, sep = "\t") %>%
        data.table::setnames("V1", "allele") %>%
        .[, type := "netMHCpan"],
      system.file("extdata",
        "mhcflurry_alleles.txt",
        package = "antigen.garnish"
      ) %>%
        data.table::fread(header = FALSE, sep = "\t") %>%
        data.table::setnames("V1", "allele") %>%
        .[, type := "mhcflurry"],
      system.file("extdata",
        "netMHCII_alleles.txt",
        package = "antigen.garnish"
      ) %>%
        data.table::fread(header = FALSE, sep = "\t") %>%
        data.table::setnames("V1", "allele") %>%
        .[, type := "netMHCII"],
      system.file("extdata",
        "netMHCIIpan_alleles.txt",
        package = "antigen.garnish"
      ) %>%
        data.table::fread(header = FALSE, sep = "\t") %>%
        data.table::setnames("V1", "allele") %>%
        .[, type := "netMHCIIpan"]
    )
  )

  # generate input for mhcflurry predictions

  mf_dt <- dt[MHC %chin% alleles[type == "mhcflurry", allele] &
    nmer_l < 15,
  .SD,
  .SDcols = c("MHC", "nmer")
  ] %>%
    data.table::copy() %>%
    data.table::setnames(c("MHC", "nmer"), c("allele", "peptide")) %>%
    unique()

  for (i in mf_dt[, allele %>% unique()]) {
    mf_dt[allele == i] %>%
      data.table::fwrite(paste0(
        "mhcflurry_input_", uuid::UUIDgenerate() %>% substr(1, 18), ".csv"
      ),
      quote = FALSE, sep = ","
      )
  }

  # generate matchable MHC substring for netMHC tools
  dt[, netMHCpan := MHC %>% stringr::str_replace(stringr::fixed("*"), "")]
  dt[, netMHC := netMHCpan %>% stringr::str_replace(stringr::fixed(":"), "")]
  dt[, netMHCII := netMHCpan %>% stringr::str_replace(stringr::fixed(":"), "") %>%
    stringr::str_replace(stringr::fixed("HLA-"), "")]
  dt[, netMHCIIpan := netMHCII %>% stringr::str_replace("DRB1", "DRB1_")]

  # replace substring with netMHC allele type
  dt[, netMHCpan := detect_mhc(netMHCpan, alleles)]
  dt[, netMHC := detect_mhc(netMHC, alleles)]
  dt[, netMHCII := detect_mhc(netMHCII, alleles)]
  dt[, netMHCIIpan := detect_mhc(netMHCIIpan, alleles)]

  dtfn <- {
    data.table::rbindlist(
      list(
        dt[!netMHC %>% is.na() &
          nmer_l < 15,
        .SD,
        .SDcols = c("netMHC", "nmer", "nmer_l")
        ] %>%
          data.table::copy() %>%
          data.table::setkey(netMHC, nmer_l) %>%
          write_netmhc_nmers("netMHC"),

        dt[!netMHCpan %>% is.na() &
          nmer_l < 15,
        .SD,
        .SDcols = c("netMHCpan", "nmer", "nmer_l")
        ] %>%
          data.table::copy() %>%
          data.table::setkey(netMHCpan, nmer_l) %>%
          write_netmhc_nmers("netMHCpan"),

        dt[!netMHCII %>% is.na() &
          nmer_l == 15,
        .SD,
        .SDcols = c("netMHCII", "nmer", "nmer_l")
        ] %>%
          data.table::copy() %>%
          data.table::setkey(netMHCII, nmer_l) %>%
          write_netmhc_nmers("netMHCII"),

        dt[!netMHCIIpan %>% is.na() &
          nmer_l == 15,
        .SD,
        .SDcols = c("netMHCIIpan", "nmer", "nmer_l")
        ] %>%
          data.table::copy() %>%
          data.table::setkey(netMHCIIpan, nmer_l) %>%
          write_netmhc_nmers("netMHCIIpan")
      )
    )
  }

  # generate commands

  dtfn[, command :=
    paste(
      type,
      "-p",
      "-l", nmer_l,
      "-a", allele,
      "-f", filename
    )]
  dtfn[type == "netMHCIIpan", command :=
    paste(
      type,
      "-inptype 1",
      "-length", nmer_l,
      "-a", allele,
      "-f", filename
    )]

  return(list(dt, dtfn))
}


#' Internal function to collate results from netMHC prediction
#'
#' @param esl List of outputs from netMHC.
#' @noRd
#' @md

collate_netMHC <- function(esl) {
  message("Collating netMHC output...")

  dtl <- lapply(
    esl, function(es) {
      command <- es[[1]]
      fn <- es[[2]]

      # read temp file
      es <- scan(file = fn, what = "character", sep = "\n")

      # check that files does not contain ERROR
      err <- es %>%
        stringr::str_detect("ERROR|Error|error") %>%
        any()
      if (err) {
        warning(paste0(command, " returned ERROR"))
        return(NULL)
      }

      # parse results
      # isolate table and header

      dtl <- es %exclude%
        "^\\#|----|Number|Distance|threshold|version|^$" %>%
        stringr::str_replace("^[ ]+", "")

      dtn <- dtl[1] %>%
        strsplit("[ ]+") %>%
        unlist() %exclude% "^$|^Bind|^Level"

      # fix table formatting
      dt <- dtl[2:length(dtl)] %>%
        stringr::str_replace_all(">", " ") %>%
        stringr::str_replace_all("=", " ") %>%
        stringr::str_replace_all("<", " ") %>%
        stringr::str_replace_all("(SB|WB)", "  ") %>%
        data.table::tstrsplit("[ ]+") %>%
        data.table::as.data.table()
      # apply names to data table
      dt %>% data.table::setnames(dt %>% names(), dtn)

      # append command
      dt$command <- command

      # set the program type from command
      ptype <- command %>% stringr::str_extract("net[A-Za-z]+")

      dtn <- dt %>% names()
      # make netMHC names consistent
      if (dtn %include% "[Pp]eptide" %>% length() > 0) {
        data.table::setnames(dt, dtn %include% "[Pp]eptide", "nmer")
      }

      if (dtn %include% "Aff.*nM.*" %>% length() > 0) {
        data.table::setnames(dt, dtn %include% "Aff.*nM.*", "affinity(nM)")
      }

      if (dtn %include% "^MHC$|HLA|Allele" %>% length() > 0) {
        data.table::setnames(dt, dtn %include% "^MHC$|HLA|Allele", "allele")
      }

      if (dtn %include% "Icore|iCore" %>% length() > 0) {
        data.table::setnames(dt, dtn %include% "Icore|iCore", "icore")
      }

      if ("Pos" %chin% dtn) {
        dt %>% data.table::setnames("Pos", "pos")
      }

      if ("Core" %chin% dtn) dt %>% data.table::setnames("Core", "core")

      # fix netMHCpan allele output to match input
      if (command %like% "netMHCpan") {
        dt[, allele := allele %>%
          stringr::str_replace(fixed("*"), "")]
      }

      # set unique column names based on program
      data.table::setnames(
        dt, dt %>% names() %exclude% "allele|nmer",
        paste0((dt %>% names() %exclude% "allele|nmer"), "_", ptype)
      )

      dt %>% data.table::setnames("allele", ptype)

      return(dt)
    }
  )
  return(dtl)
}

#' Internal function to output nmers for netMHC prediction to disk
#'
#' @param dt Data table of nmers.
#' @param type Character vector. Name of program to format for.
#'
#' @noRd
#' @md

write_netmhc_nmers <- function(dt, type) {
  if (dt %>% nrow() == 0) {
    return(NULL)
  }
  if (!c("nmer", "nmer_l") %chin% (dt %>% names()) %>% any()) {
    stop("dt must contain nmer and nmer_l columns")
  }

  # as originally written, this bugs when not all nmer lengths exist for all mhc types
  combs <- data.table::CJ(
    dt[, get(type)] %>% unique(),
    dt[, nmer_l] %>% unique()
  )

  dto <- lapply(1:nrow(combs), function(i) {
    dts <- dt[get(type) == combs$V1[i] & nmer_l == combs$V2[i]] %>%
      unique()

    if (nrow(dts) == 0) {
      return(NULL)
    }

    # parallelize over 300 peptide chunks
    chunks <- ((dts %>% nrow()) / 300) %>% ceiling()

    suppressWarnings(
      dto <- lapply(dts %>% split(1:chunks), function(dtw) {
        filename <- paste0(
          type, "_",
          uuid::UUIDgenerate() %>% substr(1, 18), ".csv"
        )


        # write out unique peptides for MHC type, length
        data.table::fwrite(dtw[, .(nmer)] %>% unique(),
          filename,
          col.names = FALSE,
          sep = ",",
          quote = FALSE
        )

        return(data.table::data.table(
          type = type,
          allele = combs$V1[i],
          nmer_l = combs$V2[i],
          filename = filename
        ))
      }) %>% data.table::rbindlist()
    )
    return(dto)
  }) %>% data.table::rbindlist()

  return(dto)
}


#' Parallelized function to create a space-separated character string between two integers
#'
#' @param x Integer. Starting integer.
#' @param y Integer. Ending integer.
#'
#' @noRd

get_ss_str <- function(x, y) {
  parallel::mcMap(function(x, y) {
    (x %>% as.integer()):(y %>% as.integer()) %>%
      paste(collapse = " ")
  }, x, y) %>% unlist()
}


#' Perform neoantigen prediction.
#'
#' Perform ensemble neoantigen prediction on a data table of missense mutations, insertions, or deletions using netMHC and mhcflurry.
#'
#' @param path Path to input `csv` or `tsv` file.
#' @param dt Data table. Input data table from `garnish_variants`, or a data table in the correct form (see [Github README](https://github.com/immune-health/antigen.garnish)).
#' @param binding_cutoff Numeric. Maximum consensus MHC-binding affinity that will be passed for IEDB and dissimilarity analysis. Default is 500 (nM). Note: If a peptide binds to any MHC allele in the table below this threshold, foreignness score and dissimilarity will be returned for all rows with that peptide.
#' @param counts Optional. A file path to a `csv` or `tsv` RNA count matrix. The first column must contain Ensembl transcript ids. All samples in the input table must be present in the count matrix.
#' @param min_counts Integer. The minimum number of estimated read counts for a transcript to be considered for neoantigen prediction. Default is 1.
#' @param peptide_length Numeric vector. Length(s) of peptides to create.
#' @param blast Logical. Run `BLASTp` to find wild-type peptide and known IEDB matches?
#' @param save Logical. Save a copy of garnish_affinity output to the working directory as "ag_output.txt"? Default is `TRUE`.
#' @param remove_wt Logical. Check all `nmer`s generated against wt peptidome and remove matches? Default is `TRUE`. If investigating wild-type sequences, set this to `FALSE`.
#' @return A data table of binding predictions including:
#' * **cDNA_seq**: mutant cDNA sequence
#' * **cDNA_locs**: starting index of mutant cDNA
#' * **cDNA_locl**: ending index of mutant cDNA
#' * **cDNA_type**: netMHC prediction tool output
#' * **frameshift**: frameshift variant?
#' * **coding**: wt cDNA sequence
#' * **coding_mut**: mutant cDNA sequence
#' * **pep_type**: type of peptide
#' * **pep_mut**: mutant peptide sequence
#' * **pep_wt**: wt peptide sequence
#' * **mismatch_s**: starting index of mutant peptide sequence
#' * **mismatch_l**: ending index of mutant peptide sequence
#' * **mutant_index**: index of mutant peptide
#' * **nmer**: nmer for prediction
#' * **nmer_i**: index of nmer in sliding window
#' * **_net**: netMHC prediction tool output
#' * **mhcflurry_**: mhcflurry_ prediction tool output
#' * **DAI**: Differential agretopicity index of missense and corresponding wild-type peptide. Differential agretopicty is the ratio of MHC binding afinity between mutant and corresponding normal peptide, with higher values indicating greater relative binding of the mutant peptide.
#' * **BLAST_A**: Ratio of consensus binding affinity of mutant peptide / closest single AA mismatch from blastp results. Returned only if `blast = TRUE`.
#'
#' antigen.garnish quality analysis metric results:
#' * **Ensemble_score**: average value of MHC binding affinity from all prediction tools.
#' * **foreignness_score**: Neoantigen foreignness threshold. Value of 0 to 1 indicating the TCR recognition probability, calculated by summing alignments in IEDB immunogenic peptides, with 1 indicating greater homology to immunogenic peptides.
#' * **IEDB_anno**: The best alignment from the IEDB database queried for the sample if applicable.
#' * **min_DAI**: Minimum of value of BLAST_A or DAI values, to provide the most conservative proteome-wide estimate of differential binding between input and wildtype matches.
#' * **dissimilarity**: Value of 0 to 1 indicating alignment to the self-proteome, calculated in an analogous manner to neoanigen foreignness, with 1 indicating greater dissimilarity.
#' @details
#' * see `list_mhc` for compatible MHC allele syntax, you may also use "all_human" or "all_mouse" in the MHC column to use all supported alleles
#'
#' Parallel cores used can be set via environment variable AG_THREADS (default: all available).
#' @seealso \code{\link{list_mhc}}
#' @seealso \code{\link{garnish_variants}}
#' @seealso \code{\link{garnish_antigens}}
#'
#' @references
#'
#' Richman LP, Vonderheide RH, and Rech AJ. Neoantigen dissimilarity to the self-proteome predicts immunogenicity and response to immune checkpoint blockade. Cell Systems. 2019.

#' Duan, F., Duitama, J., Seesi, S.A., Ayres, C.M., Corcelli, S.A., Pawashe, A.P., Blanchard, T., McMahon, D., Sidney, J., Sette, A., et al. Genomic and bioinformatic profiling of mutational neoepitopes reveals new rules to predict anticancer immunogenicity. J Exp Med. 2014.
#'
#' Luksza, M, Riaz, N, Makarov, V, Balachandran VP, et al. A neoepitope fitness model predicts tumour response to checkpoint blockade immunotherapy. Nature. 2017.
#' Rech AJ, Balli D, Mantero A, Ishwaran H, Nathanson KL, Stanger BZ, Vonderheide RH. Tumor immunity and survival as a function of alternative neopeptides in human cancer. Clinical Cancer Research, 2018.
#'
#' Wells DK, van Buuren MM, Dang KK, Hubbard-Lucey VM, Sheehan KCF, Campbell KM, Lamb A, Ward JP, Sidney J, Blazquez AB, Rech AJ, Zaretsky JM, Comin-Anduix B, Ng AHC, Chour W, Yu TV, Rizvi1 H, Chen JM, Manning P, Steiner GM, Doan XC, The TESLA Consortium, Merghoub T, Guinney J, Kolom A, Selinsky C, Ribas A, Hellmann MD, Hacohen N, Sette A, Heath JR, Bhardwaj N, Ramsdell F, Schreiber RD, Schumacher TN, Kvistborg P, Defranoux N. Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction. Cell. 2020.
#' @export garnish_affinity
#' @md

garnish_affinity <- function(dt = NULL,
                             path = NULL,
                             binding_cutoff = 500,
                             counts = NULL,
                             min_counts = 1,
                             peptide_length = 15:8,
                             blast = TRUE,
                             save = TRUE,
                             remove_wt = TRUE) {

  # generate tempdir info to use later, set early to prevent downstream errors when
  # assemble, generate, and predict steps are run separately.
  original_dir <- getwd()

  ndir <- paste("ag",
    uuid::UUIDgenerate() %>%
      stringr::str_replace_all("-", "") %>%
      substr(1, 18),
    sep = "_"
  )


  on.exit({
    message("Removing temporary files.")
    try(
      list.files(pattern = "(_nmer_fasta\\.fa)|(iedb_query.fa)|((netMHC|mhcflurry).*_[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}\\.csv)") %>% file.remove(),
      silent = TRUE
    )
    setwd(original_dir)
  })

  if (missing(dt) & missing(path)) stop("dt and path are missing.")
  if (!missing(dt) & !missing(path)) stop("Choose dt or path input.")

  if (missing(dt) & !missing(path)) {
    dt <- data.table::fread()
  }

  if (!"data.table" %chin% class(dt)) {
    stop("Input must be a data table.")
  }

  if (!"MHC" %chin% names(dt)) stop("Input must include MHC alleles, see ?garnish_affinity")

  # if class of MHC is a list column, it won't bug until first merge in make_BLAST_uuid, added this.
  if (class(dt[, MHC]) == "list") stop("MHC column must be a character column, not a list, unlist the column and rerun garnish_affinity.")

  # make sample_ids character to avoid downstream merge failure for column types
  dt[, sample_id := as.character(sample_id)]

  # remove double or more spaces in MHC string (will segfault netMHC)
  # DO NOT unique this, will not return same length vector
  dt[, MHC := MHC %>% stringr::str_replace_all("\\ +", " ")]

  input_type <- vector()

  # specify transcript vs. direct cDNA / mutant index input
  if (c("sample_id", "transcript_id", "cDNA_change", "MHC") %chin%
    (dt %>% names()) %>% all()) {
    input_type <- "transcript"
  }

  if (c("sample_id", "pep_mut", "MHC") %chin%
    (dt %>% names()) %>% all()) {
    input_type <- "peptide"
  }
  dt[, frameshift := FALSE]


  if (!
  (input_type == "transcript" ||
    input_type == "peptide")
  ) {
    stop("Incorrect input data format, see ?garnish_antigens")
  }

  if (input_type == "transcript") {
    message("Generating metadata.")

    dt %<>% get_metadata

    if ("chromosome_name" %chin% (dt %>% names())) {
      if (
        (
          dt[chromosome_name %>% stringr::str_detect(stringr::fixed(":MT"))] %>% nrow()
        ) > 0
      ) {
        warning("Removing MT transcripts.")
        dt %<>% .[!chromosome_name %>% stringr::str_detect(stringr::fixed(":MT"))]
      }
    }

    if (!missing(counts)) {
      ct <- data.table::fread(counts)

      col <- ct[, .SD, .SDcols = 1] %>% unlist()

      if (!all(col %>% stringr::str_detect(pattern = "ENS|NM_"))) {
        stop("Transcript expression matrix file first column must contain transcript ids.")
      }

      if (any(is.na(col)) |
        length(unique(col)) != length(col)) {
        stop("Transcript expression matrix has non-unique or NA values.")
      }

      ct %>% setnames(names(ct)[1], "transcript_id")

      ct %<>% melt(
        id.vars = "transcript_id",
        variable.name = "sample_id",
        value.name = "counts",
        variable.factor = FALSE
      )

      if (any(!dt[, sample_id %>% unique()] %chin% ct[, sample_id %>% unique()])) {
        stop("Transcript expression matrix does not contain columns for all samples in input data.")
      }

      ct[, counts := counts > min_counts]

      ct <- ct[counts == TRUE]

      dt <- merge(dt, ct[, .SD %>% unique(), .SDcols = c("sample_id", "transcript_id")],
        by = c("sample_id", "transcript_id")
      )

      if (nrow(dt) == 0) {
        stop("No variants in any sample met RNA counts matrix threshold, check RNA count matrix input and transcript ids or run without a count matrix.")
      }
    }


    # extract cDNA changes and make cDNA
    message("Extracting cDNA.")
    dt %<>% extract_cDNA
    message("Make cDNA.")
    dt %<>% make_cDNA

    # translate protein sequences
    dt[, pep_wt := coding %>% translate_cDNA()]
    dt[, pep_mut := coding_mut %>% translate_cDNA()]

    if (dt %>% nrow() == 0) {
      stop("No mutant peptides exist in the data table for affinity prediction.")
    }

    dt <- dt[!is.na(pep_wt) & !is.na(pep_mut)]

    dt[, frameshift := FALSE]

    dt[, cDNA_delta := ((coding_mut %>% nchar()) - (coding %>% nchar()))]

    # frameshifts have cDNA delta
    # not divisible by 3
    # does not handle stop codon loss
    # (this is acceptable because
    # no cDNA exists to know the readthrough)
    dt[cDNA_delta %% 3L != 0L, frameshift := TRUE]

    # remove variants with translated sequence-ensembl mismatch
    dt %<>% .[peptide == pep_wt]
    dt %<>% .[peptide == pep_wt]

    # remove stop codons
    # if first character is * then returns character(0) which has length zero unlike NA, so length not preserved with unlist in upcoming for loop
    # fix with:
    dt %<>% .[!stringr::str_detect(pep_mut, "^\\*")]

    for (i in dt %>% names() %include% "^pep") {
      set(dt, j = i, value = dt[, get(i) %>%
        stringr::str_extract_all("^[^\\*]+") %>%
        unlist()])
    }

    ## ---- create mutant peptide index

    message("Generating mutant peptide index.")
    # index first mismatch
    suppressWarnings({
      dt[, mismatch_s :=
        {
          (pep_wt %>%
            strsplit(split = "", fixed = TRUE) %>%
            unlist()) !=
            (pep_mut %>%
              strsplit(split = "", fixed = TRUE) %>%
              unlist())
        } %>%
        which() %>% .[1], by = 1:nrow(dt)]
    })

    # if pep_wt length < pep_mut ie stop lost variant, this returns NA so:
    # this only matters if  cDNA from meta table continues past stop (theoretical)
    dt[
      is.na(mismatch_s) & nchar(pep_wt) < nchar(pep_mut),
      mismatch_s := nchar(pep_wt) + 1
    ]

    # remove rows with matching peptides
    dt %<>% .[!pep_wt == pep_mut]

    # if pep_wt is longer than pep_mut (stop gained inframe) and no other sequence difference
    # pep_mut will be recycled and mismatch_s will not be NA
    # remove this with fixed pattern of pep_mut in pep_wt
    dt %<>% .[
      nchar(pep_wt) > nchar(pep_mut),
      drop := stringr::str_detect(pattern = stringr::fixed(pep_mut), pep_wt)
    ] %>%
      .[is.na(drop) | drop != TRUE] %>%
      .[, drop := NULL]

    # initialize mismatch length
    dt[, mismatch_l := mismatch_s]

    # frameshifts are mutants until STOP
    dt[
      frameshift == TRUE,
      mismatch_l := pep_mut %>% nchar()
    ]

    # create mutant register for
    # non-frameshift insertions
    dt[
      frameshift == FALSE & nchar(pep_mut) > nchar(pep_wt),
      mismatch_l := mismatch_s + (pep_mut %>% nchar()) -
        (pep_wt %>% nchar()) - 1
    ]

    # deletions and missense are mutants over site only
    dt[
      frameshift == FALSE & nchar(pep_mut) <= nchar(pep_wt),
      mismatch_l := mismatch_s
    ]

    # create a space-separated vector of mutant indices
    dt[, mutant_index := mismatch_s %>% as.character()]
    dt[mismatch_l > mismatch_s, mutant_index :=
      get_ss_str(mismatch_s, mismatch_l)]
  }

  if (input_type == "peptide") {
    message("Checking peptides.")

    if (any(dt[, !pep_mut %like% "^[ARNDCQEGHILKMFPSTWYV]+$"])) {
      stop("Non-standard AA one-letter codes detected in \"pep_mut\". Please remove these lines from input.")
    }

    if ("pep_wt" %chin% names(dt)) {
      if (any(dt[, !pep_wt %like% "^[ARNDCQEGHILKMFPSTWYV]+$"])) {
        stop("Non-standard AA one-letter codes detected in \"pep_wt\". Please remove these lines from input.")
      }

      if (any(dt[, mutant_index %>% unique()] %like% "\\ ")) {
        stop(paste(
          "MNV and frameshifts are not supported in paired mutant wild-type peptide input mode.",
          "Please provide a single amino acid position as \"mutant_index\" or use \"pep_mut\" input only.",
          sep = "\n"
        ))
      }

      if (any(dt[, stringr::str_detect(pattern = stringr::fixed(pep_mut), stringr::fixed(pep_wt))])) {
        warning(
          paste("Rows where pep_mut contained no mutant sequence have been dropped:",
            paste(
              dt[, which(stringr::str_detect(pattern = stringr::fixed(pep_mut), stringr::fixed(pep_wt)))],
              collapse = ", "
            ),
            sep = "\n"
          )
        )

        dt <- dt[!stringr::str_detect(pattern = stringr::fixed(pep_mut), stringr::fixed(pep_wt))]

        if (nrow(dt) == 0) {
          return("no variants for peptide generation")
        }
      }
    }

    dt[mutant_index == "all", mutant_index :=
      get_ss_str(1, pep_mut %>% nchar())]
  }

  message("Generating variants")

  # generation a uuid for each unique variant

  suppressWarnings(dt[, var_uuid := uuid::UUIDgenerate(n = .N)])

  # separate over mutant indices

  if (dt[, mutant_index %>% unique()] %>%
    stringr::str_detect(" ") %>% any()) {
    dts <- dt %>% tidyr::separate_rows("mutant_index", sep = " ")
  } else {
    dts <- dt
  }

  dts %<>% data.table::as.data.table(.)

  # convert back to numeric

  dts[, mutant_index := mutant_index %>% as.numeric()]

  # generate a data table of unique variants for peptide generation
  # this gives index error with frameshift inputs, only occurs when separate_rows is used above
  # open issue on data.table repo, name exists but column index is invalid:
  # Error in getindex(x, names(x)[xcols]) :
  # Internal error: index 'frameshift' exists but is invalid
  # I randomly fixed by setting key to another column, at some point key was set
  # to frameshift, this produces error on dtfs and dtnfs creation
  # data.table::copy() at dts creation did not solve this
  setkey(dts, "sample_id")

  dtnfs <- dts[frameshift == FALSE]
  dtfs <- dts[frameshift == TRUE]

  if (input_type == "transcript") {
    basepep_dt <- data.table::rbindlist(list(
      # take pep_wt for non-fs for DAI calculation
      data.table::rbindlist(list(
        dtnfs %>%
          data.table::copy() %>%
          .[, pep_base := pep_wt] %>%
          .[, pep_type := "wt"],
        dtnfs %>%
          data.table::copy() %>%
          .[, pep_base := pep_mut] %>%
          .[, pep_type := "mutnfs"]
      )),
      # take only pep_mut for fs
      dtfs %>%
        data.table::copy() %>%
        .[, pep_base := pep_mut] %>%
        .[, pep_type := "mut_other"]
    )) %>%

      # take unique peptides
      .[, .SD, .SDcols = c(
        "var_uuid",
        "pep_type",
        "pep_base",
        "mutant_index"
      )] %>%
      unique()
  }

  if (input_type == "peptide" & !"pep_wt" %chin% names(dtnfs)) {
    basepep_dt <- dtnfs %>%
      data.table::copy() %>%
      .[, pep_base := pep_mut] %>%
      .[, pep_type := "mut_other"]
  }

  if (input_type == "peptide" & "pep_wt" %chin% names(dtnfs)) {
    basepep_dt <- data.table::rbindlist(list(
      # take pep_wt for non-fs for DAI calculation
      dtnfs %>%
        data.table::copy() %>%
        .[, pep_base := pep_wt] %>%
        .[, pep_type := "wt"],
      dtnfs %>%
        data.table::copy() %>%
        .[, pep_base := pep_mut] %>%
        .[, pep_type := "mutnfs"]
    ))
  }

  # filter mutant_index
  basepep_dt <- basepep_dt[(nchar(pep_base) - mutant_index) >= 0]

  if (basepep_dt %>% nrow() == 0) {
    return("no variants for peptide generation")
  }

  nmer_dt <- make_nmers(basepep_dt, peptide_length)

  nmer_dt %<>% .[, nmer_l := nmer %>% nchar()]

  dt <- merge(dt, nmer_dt,
    by = intersect(names(dt), names(nmer_dt)),
    all.x = TRUE
  )

  # remove peptides present in global normal protein database

  # load global peptide database
  if (remove_wt) {
    message("Filtering WT peptide matches.")

    d <- system.file(package = "antigen.garnish") %>% file.path(., "extdata")

    if (
      (!dt$MHC %like% "HLA" %>% any()) &
        (!dt$MHC %like% "H-2" %>% any()) &
        !all(dt$MHC %chin% c("all_human", "all_mouse"))
    ) {
      stop("MHC do not contain \"HLA-\" or \"H-2\" as a pattern.
              Alleles must be correctly formatted, see list_mhc().")
    }

    if (
      (dt$MHC %like% "HLA" %>% any() &
        !dt$MHC %like% "H-2" %>% any()) ||
        (dt$MHC == "all_human") %>% any()
    ) {
      pepv <-
        file.path(
          d,
          "antigen.garnish_GRCh38_pep.RDS"
        ) %>%
        readRDS(.)
    }
    if (
      (dt$MHC %like% "H-2" %>% any() &
        !dt$MHC %like% "HLA" %>% any()) ||
        (dt$MHC == "all_mouse") %>% any()
    ) {
      pepv <-
        file.path(
          d,
          "antigen.garnish_GRCm38_pep.RDS"
        ) %>%
        readRDS(.)
    }

    if (
      (dt$MHC %like% "HLA" %>% any() &
        dt$MHC %like% "H-2" %>% any()) ||
        (any(dt$MHC == "all_human") & any(dt$MHC == "all_mouse"))
    ) {
      pepv <- c(
        file.path(d, "antigen.garnish_GRCh38_pep.RDS") %>%
          readRDS(.),
        file.path(d, "antigen.garnish_GRCm38_pep.RDS") %>%
          readRDS(.)
      )
    }

    # get unique normal proteins and unique non-wt nmers to match

    pepv %<>% unique

    nmv <- dt[pep_type != "wt", nmer %>% unique()]

    # return vector of matched nmers to drop

    mv <- lapply(nmv %>% seq_along(), function(x) {
      ifelse(
        stringr::str_detect(pepv, stringr::fixed(nmv[x])) %>% any(),
        return(nmv[x]),
        return(NULL)
      )
    }) %>% unlist()

    # drop matched nmers
    dt %<>% .[!(nmer %chin% mv & pep_type != "wt")]
  }

  # generation a uuid for each unique nmer

  nmer_uuid_dt <- dt[, .SD, .SDcols = "nmer"] %>% unique()

  nmer_uuid_dt[, nmer_uuid :=
    uuid::UUIDgenerate(use.time = FALSE, n = .N)]

  dt %<>% merge(nmer_uuid_dt, by = "nmer", all.x = TRUE)


  if (input_type == "transcript" ||
    (input_type == "peptide" & "pep_wt" %chin% names(dt))
  ) {
    dt %<>% make_DAI_uuid

    # remove mut == wt by sample_id
    for (id in dt[, sample_id %>% unique()]) {

      # get wt nmers for fast !chmatch
      id_wt_nmers <- dt[
        sample_id == id &
          pep_type == "wt",
        nmer %>% unique()
      ]

      dt %<>% .[pep_type == "wt" |
        sample_id != id |
        (sample_id == id &
          pep_type != "wt" &
          !nmer %chin% id_wt_nmers)]
    }
  }

  # DAI cannot be calculated with peptide input
  if (input_type == "peptide" & !"pep_wt" %chin% names(dt)) {
    dt[, dai_uuid := NA %>% as.character()]
  }

  if (blast) {
    dt %<>% make_BLAST_uuid
  }

  dir.create(ndir)

  setwd(ndir)

  message("Generating prediction commands.")
  dtl <- dt %>% get_pred_commands()

  # run commands

  if (check_pred_tools() %>% unlist() %>% all()) {
    dto <- run_netMHC(dtl[[2]])
    run_mhcflurry()
    dt <- merge_predictions(dto, dtl[[1]])
  } else {
    warning("Missing prediction tools in PATH, returning without predictions.")
  }

  message("Setting up dissimilarity calculation.")

  # now  that we know binders, get our foreignness_score and dissimilarity values
  # iterate over multiple species
  dt[MHC %like% "HLA", spc := "human"]
  dt[MHC %like% "H-2", spc := "mouse"]

  mns <- dt[pep_type != "wt" & Ensemble_score < binding_cutoff & spc == "mouse", nmer %>% unique()]
  hns <- dt[pep_type != "wt" & Ensemble_score < binding_cutoff & spc == "human", nmer %>% unique()]

  if (length(mns) != 0 & blast) {
    idt <- mns %>%
      foreignness_score(db = "mouse") %>%
      .[, spc := "mouse"]

    sdt <- mns %>%
      dissimilarity_score(db = "mouse") %>%
      .[, spc := "mouse"]

    dt <- merge(dt, idt, all.x = TRUE, by = c("nmer", "spc"))

    dt <- merge(dt, sdt, all.x = TRUE, by = c("nmer", "spc"))
  }

  if (length(hns) != 0 & blast) {
    idt <- hns %>%
      foreignness_score(db = "human") %>%
      .[, spc := "human"]

    sdt <- hns %>%
      dissimilarity_score(db = "human") %>%
      .[, spc := "human"]

    dt <- merge(dt, idt, all.x = TRUE, by = c("nmer", "spc"))

    dt <- merge(dt, sdt, all.x = TRUE, by = c("nmer", "spc"))
  }

  dt[, spc := NULL]

  # set NA to 0 for foreignness_score to match original implementation
  # set NA to 0 for dissimilarity, no alignments is truly dissimilar
  if ("foreignness_score" %chin% names(dt)) {
    dt[is.na(foreignness_score) & pep_type != "wt" & Ensemble_score < binding_cutoff, foreignness_score := 0]
  }

  if ("dissimilarity" %chin% names(dt)) {
    dt[is.na(dissimilarity) & pep_type != "wt" & Ensemble_score < binding_cutoff, dissimilarity := 0]
  }

  if (blast) {

    # calculate minimum proteome-wide DAI
    if ("blast_uuid" %chin% names(dt)) {
      dt[!is.na(BLAST_A), min_DAI := BLAST_A]
    }

    if ("DAI" %chin% names(dt)) {

      # edge case with small table or deliberately not running blast
      if (!"min_DAI" %chin% names(dt)) {
        dt[, min_DAI := DAI]
      }

      v <- 1:nrow(dt[!is.na(DAI)])

      if (nrow(dt[!is.na(DAI)]) != 0) {
        dt[!is.na(DAI), min_DAI := min(DAI, min_DAI, na.rm = TRUE),
          by = v
        ]
      }
    }
  }

  gplot_fn <- format(Sys.time(), "%d/%m/%y %H:%M:%OS") %>%
    stringr::str_replace_all("[^A-Za-z0-9]", "_") %>%
    stringr::str_replace_all("[_]+", "_")

  if (save) {
    fn <- paste0("ag_output_", gplot_fn, ".txt")
    dt %>% data.table::fwrite(fn,
      sep = ",",
      quote = FALSE,
      row.names = FALSE
    )

    file.copy(fn, file.path(original_dir, fn))

    if (file.exists(file.path(original_dir, fn))) {
      setwd(original_dir)
      if (dir.exists(ndir)) {
        unlink(ndir, recursive = TRUE, force = TRUE)
      }
    }
  }

  if (!save) {
    setwd(original_dir)

    if (dir.exists(ndir)) {
      unlink(ndir, recursive = TRUE, force = FALSE)
    }
  }

  return(dt)
}


#' Internal function for parallelized `nmer` creation.
#'
#' @param dt Data table. Input data table from `garnish_affinity`.
#' @param plen Numeric vector. Length(s) of peptides to create.
#'
#' @noRd
#' @md

make_nmers <- function(dt, plen = 8:15) {
  if (!plen %>% is.numeric()) {
    stop("peptide length must be numeric")
  }

  if (!(plen %in% 8:15 %>% any())) {
    stop("peptide length must be numeric")
  }

  if (!c(
    "var_uuid",
    "pep_type",
    "pep_base",
    "mutant_index"
  ) %chin% (dt %>% names()) %>% all()) {
    stop("dt is missing columns")
  }

  dt %<>% data.table::as.data.table()

  lines <- nrow(dt)

  # for every peptide length from greatest to least
  plen %<>% sort %>% rev()


  # a function to generate sliding
  # window n-mers over an index of interest

  message("Generating nmers")


  .fn <- purrr::partial(paste, collapse = "")

  nmer_list <- lapply(
    1:nrow(dt),

    function(n) {

      ## --- Write peptide fragments

      nmer_list <- lapply(plen, function(pl) {
        mut_frag_t <- unlist(strsplit(dt$pep_base[n], "", fixed = TRUE))
        mut_frag_loc <- dt$mutant_index[n]

        # if the peptide is not long enough, return
        if (!(length(mut_frag_t)) >= pl) {
          return(NULL)
        }

        # re-register peptide if the mutant index
        # is not centered due to back truncation

        # trim peptide front to mutant site
        if (mut_frag_loc > pl) {
          rml <- mut_frag_loc - pl + 1
          mut_frag_t <- mut_frag_t[rml:(length(mut_frag_t))]
          mut_frag_loc <- pl
        }

        # trim peptide back to mutant site
        fpl <- mut_frag_loc + pl - 1
        if (fpl < length(mut_frag_t)) {
          mut_frag_t <- mut_frag_t[1:(fpl)]
        }

        # slide across the peptide window
        # create (n = pl )-mers wrapped in
        # sync to prevent zoo::rollapply() stdout

        nmers <- zoo::rollapply(mut_frag_t, pl,
          .fn,
          partial = FALSE,
          align = "left"
        )


        # return a data table of peptides

        return(
          data.table::data.table(
            nmer = nmers,
            nmer_i = 1:length(nmers),
            pep_base = dt$pep_base[n],
            pep_type = dt$pep_type[n],
            var_uuid = dt$var_uuid[n]
          )
        )
      })

      nmer_dt <- data.table::rbindlist(nmer_list)

      return(nmer_dt)
    }
  )

  nmer_dt <- data.table::rbindlist(nmer_list)
  return(nmer_dt)
}
