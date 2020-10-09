#' Internal function to create Ensembl transcript database
#'
#' Internal function to create the antigen.garnish data directory Ensembl portion of transcript database file `GRChm38_meta.RDS` from [GRCh38.cds.all.fa](ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/) Ensembl files.
#'
#'
#' @noRd
#' @md

make_ensembl_db <- function() {

  # load remaining old Hu and Ms databases and add transcripts if they do not exist in newer databases
  # thus, we prioritize the current database
  ret <- list.files(pattern = "GRC.*38.*fa") %>%
    lapply(function(fn) {
      dt <- data.table::fread(fn, header = FALSE, sep = "\t")

      # index of FASTA entries
      nameIdx <- dt[, (V1 %likef% "ENS") %>% which()] # header lines
      # add final index for EOF
      nameIdx %<>% c(nrow(dt) + 1)

      ret <- 1:(length(nameIdx) - 1) %>% # every row except last, which is EOF + 1
        parallel::mclapply(function(i) {
          data.table::data.table(
            description = dt[nameIdx[i]], # header is the first line
            coding = dt[(
              # sequence is to next index (or EOF - 1 = EOF)
              # start: index + 1
              # end: next start - 1
              (nameIdx[i] + 1):(nameIdx[i + 1] - 1)
            )]$V1 %>%
              paste(collapse = "")
          )
        }) %>% data.table::rbindlist(use.names = TRUE, fill = TRUE)
    }) %>%
    data.table::rbindlist(.)

  ret %>% data.table::setnames("description.V1", "description")

  # parse additional columns
  ret[, chromosome_name := description %>% stringr::str_extract("(?<=chromosome:)[^ ]+")]
  ret[, ensembl_gene_id := description %>% stringr::str_extract("(ENSMUSG|ENSG)[0-9]+(\\.[0-9]+)?")]
  ret[, transcript_id := description %>% stringr::str_extract("(ENSMUST|ENST)[0-9]+(\\.[0-9]+)?")]

  # take only unique older transcripts
  ret %<>% unique(by = "transcript_id")

  # do not use non-versioned transcripts because they are not reliable
  ret[transcript_id %likef% "."] %>% nrow()
  ret[!transcript_id %likef% "."] %>% nrow()
  ret %<>% .[transcript_id %likef% "."]

  # check integrity of parsing
  ret[!ensembl_gene_id %likef% "ENS"] %>% nrow() == 0
  ret[!transcript_id %likef% "ENS"] %>% nrow() == 0
  ret[coding %likef% "ENS"] %>% nrow() == 0

  # remove MT
  ret %<>% .[!chromosome_name %likef% ":MT"]

  # predict peptides from mRNA
  uniqCoding <- ret[, .SD, .SDcols = "coding"] %>% unique()
  uniqCoding[, peptide := coding %>% parallel::mclapply(function(i) {
    i %>%
      Biostrings::DNAString() %>%
      Biostrings::translate(no.init.codon = TRUE, if.fuzzy.codon = "solve") %>%
      as.character(.)
  }) %>% unlist()]

  final <- merge(ret, uniqCoding, by = "coding", all.x = TRUE)

  final %>% saveRDS("GRChm38_meta.RDS", compress = FALSE)
}


#' Internal function to create RefSeq transcript database
#'
#' Internal function to create the antigen.garnish data directory Refseq portion of transcript database file `GRChm38_meta.RDS`.
#'
#' @param txt Text file with one Refseq transcript ID (e.g. NM_000535.5) per row, compiled from current RefSeq database pull.
#'
#' @return A data table of Refseq transcript IDs with metadata.
#'
#' @noRd
#' @md

make_refseq_db <- function(txt) {
  rs <- data.table::fread(txt)

  rs[, refseq_id_long := raw %>% stringr::str_extract("NM_[0-9]+\\.[0-9]+")]
  rs[, refseq_id := refseq_id_long %>% stringr::str_extract("NM_[0-9]+")]
  rs[, version := refseq_id_long %>% stringr::str_replace("NM_[0-9]+\\.", "") %>% as.numeric()]

  # for each versioned identifier, create string
  # of all previous identifiers
  rs[, refseq_id_long := paste0(refseq_id, ".", (version - 1):1) %>%
    paste(collapse = " "), by = 1:nrow(rs)]

  # separate each identifer to a row
  rsAllIds <- tidyr::separate_rows(rs, refseq_id_long, sep = " ") %>%
    data.table::as.data.table() %>%
    unique(by = "refseq_id_long")

  # save output containing all possible refseq versions
  rsAllIds %>% data.table::fwrite("20200914-refseq_mrna_versions.txt", sep = "\t")

  out <- rsAllIds[, refseq_id_long]

  l <- out %>% split(1:500)



  # query NIH eutils via http for transcript annotations
  1:length(l) %>% lapply(function(i) {
    print(i)
    ids <- l[[i]] %>% paste(collapse = ",")
    glue::glue("http 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ids}&retmode=xml' > rs{i}.xml") %>% system()
  })

  # extract RefSeq IDm cDNA,  and peptide sequence from .xml
  dt <- list.files(pattern = "^[0-9]+\\.xml$") %>%
    parallel::mclapply(function(fN) {
      print(fN)
      x <- xml2::read_xml(fN)

      # extract identifer and protein translation
      refseq_id_long <- x %>%
        xml2::xml_find_all("//GBInterval_accession") %>%
        xmltools::xml_dig_df(dig = TRUE) %>%
        unlist() %>%
        unique() %>% # listed multiple times in each XML entry
        data.table::as.data.table(.) %>%
        data.table::setnames("refseq_id_long")

      ret <- refseq_id_long[, refseq_id_long] %>%
        lapply(function(id) {
          if (!id %like% "NM") {
            return(NULL)
          }

          sub <- x %>%
            # move to sequence record ancestor corresponding to transcript ID
            xml2::xml_find_all(glue::glue("//GBSeq_accession-version[text()='{id}']/ancestor::GBSeq"))

          peptide <- sub %>%
            # find translation entry
            xml2::xml_find_all(".//GBQualifier_name[text()='translation']") %>%
            # extract first and only sibling that contains peptide sequence
            xml2::xml_find_all("./following-sibling::*") %>%
            xml2::xml_text(.)

          codingTo <- sub %>%
            # find translation entry
            xml2::xml_find_all(".//GBFeature_key[text()='CDS']") %>%
            # extract first and only sibling that contains peptide sequence
            xml2::xml_find_all("./following-sibling::*") %>%
            xml2::xml_find_all(".//GBInterval_to") %>%
            xml2::xml_text(.) %>%
            as.numeric()

          codingFrom <- sub %>%
            # find translation entry
            xml2::xml_find_all(".//GBFeature_key[text()='CDS']") %>%
            # extract first and only sibling that contains peptide sequence
            xml2::xml_find_all("./following-sibling::*") %>%
            xml2::xml_find_all(".//GBInterval_from") %>%
            xml2::xml_text(.) %>%
            as.numeric()

          cDNA <- sub %>%
            # find translation entry
            xml2::xml_find_all(".//GBSeq_sequence") %>%
            xml2::xml_text(.) %>%
            toupper()

          # return NA if transcript information is not available
          # otherwise, return transcript information
          if (
            (c(
              length(peptide),
              length(codingTo),
              length(codingFrom),
              length(cDNA)
            ) == 0) %>% any()
          ) {
            ret <- data.table(
              refseq_id_long = id,
              peptide = as.character(NA),
              codingTo = as.character(NA),
              codingFrom = as.character(NA),
              cDNA = as.character(NA),
              coding = as.character(NA),
              pepFromCoding = as.character(NA)
            )
          } else {
            coding <- cDNA %>% substr(codingFrom, codingTo) # get coding portion of cDNA
            pepFromCoding <- coding %>%
              translate_cDNA(.) %>%
              stringr::str_replace("\\*", "")
            ret <- data.table(
              refseq_id_long = id, # ID
              peptide = peptide, # reference peptide
              codingTo = codingTo, # cDNA coding start position
              codingFrom = codingFrom, # cDNA coding end position
              cDNA = cDNA, # reference cDNA sequence
              coding = coding, # CDS squence
              pepFromCoding = pepFromCoding # test antigen.garnish peptide translation
            )
          }

          return(ret)
        }) %>%
        data.table::rbindlist(.)
    }) %>%
    data.table::rbindlist(.)

  return(dt)
}

#' Internal function to determine if a RefSeq CDS has changed across transcript versions
#'
#' @param id Character vector of Refseq transcript IDs (e.g. NM_000535.5)
#'
#' @return A data table of Refseq transcript IDs with logical columns indicating if the CDS matches and previous version of the transcript or the next version of the transcript.
#'
#' @noRd
#' @md

checkRefSeqVersion <- function(ids) {
  ag_dirs <- c(
    paste0(
      getwd(),
      "/",
      "antigen.garnish"
    ),
    paste0(
      Sys.getenv("HOME"),
      "/",
      "antigen.garnish"
    )
  )

  if (Sys.getenv("AG_DATA_DIR") == "") {
    for (i in ag_dirs) {
      if (i %>% dir.exists()) {
        ag_dir <- i
        break
      }
    }
  }

  db <- file.path(ag_dir, "/GRChm38_meta.RDS") %>%
    readRDS(.)

  ret <- ids %>%
    lapply(function(id) {
      version <- id %>%
        stringr::str_replace("NM_[0-9]+\\.", "") %>%
        as.numeric()

      if (version %>% is.na() || length(version) == 0) {
        stop("cannot determine transcript version")
      }

      oldId <- paste0(
        id %>% stringr::str_extract("NM_[0-9]+"), ".",
        version - 1
      )

      newId <- paste0(
        id %>% stringr::str_extract("NM_[0-9]+"), ".",
        version + 1
      )

      browser()

      codingOld <- db[transcript_id %chin% c(id, oldId), coding]
      codingNew <- db[transcript_id %chin% c(id, newId), coding]

      dt <- data.table::data.table(
        transcript_ida = id,
        previousVersionCDSMatches = codingOld[1] == codingOld[2],
        nextVersionCDSMatches = codingNew[1] == codingNew[2]
      )
      return(dt)
    }) %>%
    data.table::rbindlist(.)

  return(ret)
}
