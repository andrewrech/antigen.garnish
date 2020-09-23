#' Internal function to create Ensembl transcript dataabse
#'
#' Internal function to create the antigen.garnish data directory Ensembl transcript database file `GRChm38_meta.RDS` from [GRCh38.cds.all.fa](ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/) Ensembl files.
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
