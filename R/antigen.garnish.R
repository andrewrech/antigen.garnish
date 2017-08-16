#' \pkg{antigen.garnish}: ensemble neoepitope prediction from DNA variants in R.
#'
#' [antigen.garnish](http://neoepitopes.io) is an R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA point mutations, insertions, and deletions in VCF format and performs neoepitope prediction. Output is individual peptides and a summary by sample.
#'
#'Advantages
#'
#'1. **Simplicity**: summarized neoepitopes for each sample
#'1. **Thoroughness**:
#'    - missense mutations and frameshifts
#'    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i).
#'1. **Speed**:
#'    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
#'    - on an Amazon Web Services `m4.16xlarge` EC2 instance, 20,000 consensus predictions using 100+ MHC types in under 5 minutes
#' @section Manifest:
#' * `garnish_variants`: process variants from [SnpEff](http://snpeff.sourceforge.net/)
#' * `garnish_predictions`: perform ensemble neoepitope prediction
#' * `garnish_summary`: summarize neoepitope prediction
#' @docType package
#' @name antigen.garnish
#' @import colorspace
#' @import data.table
#' @import dt.inflix
#' @import parallel
#' @import stringr
#' @import testthat
#' @importFrom Biostrings DNAString translate
#' @importFrom biomaRt useMart getBM getSequence
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom stats na.omit
#' @importFrom stringi stri_detect_fixed
#' @importFrom tidyr separate_rows
#' @importFrom utils download.file
#' @importFrom uuid UUIDgenerate
#' @importFrom vcfR read.vcfR
#' @importFrom zoo rollapply
#' @md
NULL



## ---- get_DAI_uuid
#' Internal function to pair peptides from missense sites by a common UUID for DAI calculations
#'
#' @param dt Data table of nmers.
#'
#' @export get_DAI_uuid
#' @md
get_DAI_uuid <- function(dt){

  if (!c("pep_type", "nmer", "nmer_i",
          "nmer_l", "var_uuid", "frameshift",
          "pep_mut", "pep_wt") %chin%
     (dt %>% names) %>% any)
     stop("dt is missing columns")

    daidt <- dt %>%
        # explicitly select missense
        .[frameshift == FALSE &
        (pep_mut %>% nchar) == (pep_wt %>% nchar)] %>%
     # recombine table, creating pairs
     # merge vs. sort and bind so edges cases
        # result in lost peptides vs. total mis-
        # alignment
      {
     merge(
       .[pep_type == "wt", .SD,
       .SDcols = c("nmer", "nmer_i",
                   "nmer_l", "var_uuid")] %>%
       data.table::setnames("nmer", "wt_nmer"),

       .[pep_type == "mutnfs", .SD,
       .SDcols = c("nmer", "nmer_i",
                   "nmer_l", "var_uuid")] %>%
       data.table::setnames("nmer", "mtnfs_nmer"),
     by = c("var_uuid", "nmer_i", "nmer_l"))
      } %>%
      # create a DAI uuid
     .[, dai_uuid := parallel::mclapply(1:nrow(.),
                     uuid::UUIDgenerate) %>% unlist] %>%

      # bind back into one table
      {
        rbindlist(list(
              .[, .SD, .SDcols = c("wt_nmer", "nmer_i",
                 "nmer_l", "var_uuid", "dai_uuid")] %>%
              data.table::setnames("wt_nmer", "nmer"),
              .[, .SD, .SDcols = c("mtnfs_nmer", "nmer_i",
                 "nmer_l", "var_uuid", "dai_uuid")] %>%
              data.table::setnames("mtnfs_nmer", "nmer")
                     ))
      }
      # merge back together
      dt %<>% merge(daidt,
                    by = c("nmer",
                           "nmer_i",
                           "nmer_l",
                           "var_uuid"),
                    all.x = TRUE)
      return(dt)
}



## ---- merge_predictions
#' Internal function to merge input data table and prediction results from netMHC tools, mhcflurry
#'
#' @param l Output list from run_netMHC
#' @param dt Input data table.
#' @export merge_predictions
#' @md
 merge_predictions <- function(l, dt){

      message("Merging output")
      # merge netMHC by program type
        progl <- lapply(l %>% seq_along, function(dti) {
         l[[dti]]$command[1] %>% stringr::str_extract("net[A-Za-z]+")
         })

       # merge netMHC output

        for (ptype in (progl %>% unique %>% unlist)) {

          dt <- merge(dt, l[(progl == ptype) %>% which] %>%
                data.table::rbindlist, by = c("nmer", ptype), all.x = TRUE,
                allow.cartesian = TRUE)
        }

      # read and merge mhcflurry output
      fdt <- list.files(pattern = "mhcflurry_output.*csv") %>%
        lapply(., function(x){
          data.table::fread(x)
        }) %>% data.table::rbindlist %>%
            data.table::setnames(c("allele", "peptide"), c("MHC", "nmer"))
            dt <- merge(dt, fdt, by = c("nmer", "MHC"), all.x = TRUE)

      # calculate netMHC consensus score, preferring non-*net tools
          for (col in (dt %>% names %include% "aff|[Rr]ank|Consensus_scores")) {
            suppressWarnings(set(dt, j = col, value = dt[, get(col) %>% as.numeric]))
          }

      # get vector of netMHC scores
         cols <- dt %>% names %includef% c("affinity(nM)")

      # create a long format table
        # to calculate consensus score
      dtm <- dt[, .SD, .SDcols = c("nmer", "MHC", cols)] %>%
           melt(id.vars = c("nmer", "MHC")) %>%
      # order affinity predictions by program preference
           .[, variable := variable %>% factor(levels = cols)] %>%
      # key table so first non-NA value is the preferred programe
           data.table::setkey(nmer, MHC, variable) %>%
           .[, .(Consensus_scores =
      # define Consensus_score
              value %>%
              na.omit %>%
              .[1]), by = c("nmer", "MHC")] %>%
           .[!Consensus_scores %>% is.na]

      # merge back
      dt %<>% merge(dtm, by = c("nmer", "MHC"))

      # take average of mhcflurry and best available netMHC tool
      dt[, Consensus_scores := c(Consensus_scores,
                                       mhcflurry_prediction) %>%
                                      mean(na.rm = TRUE),
                                      by = 1:nrow(dt)]

        dt[, DAI := NA %>% as.numeric]

        data.table::setkey(dt, pep_type, dai_uuid)

        # dai_uuid is always length 2
        # calculate DAI by dividing

        dt[!dai_uuid %>% is.na,
        DAI := Consensus_scores[2] /
          Consensus_scores[1], by = dai_uuid]

      return(dt)
    }



## ---- get_pred_commands
#' Internal function to create commands for neoepitope prediction.
#'
#' @param dt Data.table of predictions to run.
#' @export get_pred_commands
#' @md
get_pred_commands <- function(dt){

  if (!c("nmer", "MHC", "nmer_l") %chin%
     (dt %>% names) %>% any)
     stop("dt is missing columns")

  if (dt[, MHC %>% unique] %>%
      stringr::str_detect(" ") %>% any) dt %<>%
      tidyr::separate_rows("MHC", sep = " ")

  # get available MHC alleles for predictions

  alleles <- data.table::rbindlist(
                         list(
  system.file("extdata",
      "netMHC_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHC"],
  system.file("extdata",
      "netMHCpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCpan"],
  system.file("extdata",
      "mhcflurry_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcflurry"],
  system.file("extdata",
      "netMHCII_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCII"],
  system.file("extdata",
      "netMHCIIpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCIIpan"]
                      ))

  # generate input for mhcflurry predictions

    mf_dt <- dt[MHC %chin% alleles[type == "mhcflurry", allele] &
    nmer_l < 15,
          .SD, .SDcols = c("MHC", "nmer")] %>%
    data.table::copy %>%
      data.table::setnames(c("MHC", "nmer"), c("allele", "peptide")) %>%
      unique

    for (i in mf_dt[, allele %>% unique]){
      mf_dt[allele == i] %>%
      data.table::fwrite(paste0(
            "mhcflurry_input_",  uuid::UUIDgenerate(), ".csv"))
    }

  # generate matchable MHC substring for netMHC tools
    dt[, netMHCpan := MHC %>% stringr::str_replace(stringr::fixed("*"), "")]
    dt[, netMHC := netMHCpan %>% stringr::str_replace(stringr::fixed(":"), "")]
    dt[, netMHCII := netMHCpan %>% stringr::str_replace(stringr::fixed(":"), "") %>%
                                       stringr::str_replace(fixed("HLA-"), "")]
    dt[, netMHCIIpan := netMHCII %>% stringr::str_replace("DRB1", "DRB1_")]

  # replace substring with netMHC allele type

    dt[, netMHCpan := detect_hla(netMHCpan, alleles)]
    dt[, netMHC := detect_hla(netMHC, alleles)]
    dt[, netMHCII := detect_hla(netMHCII, alleles)]
    dt[, netMHCIIpan := detect_hla(netMHCIIpan, alleles)]

    dtfn <-
      {
        data.table::rbindlist(
          list(
            dt[!netMHC %>% is.na &
            nmer_l < 15,
              .SD, .SDcols = c("netMHC", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHC, nmer_l) %>%
            write_nmers("netMHC"),

            dt[!netMHCpan %>% is.na &
            nmer_l < 15,
              .SD, .SDcols = c("netMHCpan", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHCpan, nmer_l) %>%
            write_nmers("netMHCpan"),

            dt[!netMHCII %>% is.na &
            nmer_l == 15,
              .SD, .SDcols = c("netMHCII", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHCII, nmer_l) %>%
            write_nmers("netMHCII"),

            dt[!netMHCIIpan %>% is.na &
            nmer_l == 15,
              .SD, .SDcols = c("netMHCIIpan", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHCIIpan, nmer_l) %>%
            write_nmers("netMHCIIpan")
                                   ))
      }

  # generate commands

    dtfn[, command :=
      paste(
            type,
            "-p",
            "-l", nmer_l,
            "-a", allele,
            "-f", filename)
        ]
    dtfn[type == "netMHCIIpan", command :=
      paste(
            type,
            "-inptype 1",
            "-length", nmer_l,
            "-a", allele,
            "-f", filename)
        ]
    return(list(dt, dtfn))
  }



## ---- collate_netMHC
#' Internal function to collate results from netMHC prediction
#'
#' @param esl List of outputs from netMHC.
#' @export collate_netMHC
#' @md
collate_netMHC <- function(esl){

   dtl <- parallel::mclapply(
         esl, function(es){

          command <- es[[1]]
          es <- es[[2]]

          # parse results
            # isolate table and header

              dtl <- es %exclude%
                "^\\#|----|Number|Distance|threshold|version|^$" %>%
                stringr::str_replace("^[ ]+", "")

              dtn <- dtl[1] %>%
                strsplit("[ ]+") %>%
                unlist %exclude% "^$|^Bind|^Level"

          # if error, warn and return
          if (dtn %>%
              stringr::str_detect("ERROR") %>%
              any) {
            warning(paste0(command, " returned ERROR"))
            return(NULL)
          }

           # fix table formatting
              dt <- dtl[2:length(dtl)] %>%
                stringr::str_replace_all(">", " ") %>%
                stringr::str_replace_all("=", " ") %>%
                stringr::str_replace_all("<", " ") %>%
                stringr::str_replace_all("(SB|WB)", "  ") %>%
                data.table::tstrsplit("[ ]+") %>%
                data.table::as.data.table
          # apply names to data table
            dt %>% data.table::setnames(dt %>% names, dtn)

          # append command
            dt$command <- command

          # set the program type from command
            ptype <- command %>% stringr::str_extract("net[A-Za-z]+")

          dtn <- dt %>% names
          # make netMHC names consistent
            if (dtn %include% "[Pp]eptide" %>% length > 0)
                data.table::setnames(dt, dtn %include% "[Pp]eptide", "nmer")

            if (dtn %include% "Aff.*nM.*" %>% length > 0)
                data.table::setnames(dt, dtn %include% "Aff.*nM.*", "affinity(nM)")

            if (dtn %include% "HLA|Allele" %>% length > 0)
                data.table::setnames(dt, dtn %include% "HLA|Allele", "allele")

            if (dtn %include% "Icore|iCore" %>% length > 0)
                data.table::setnames(dt,dtn %include% "Icore|iCore", "icore")

            if ("Pos" %chin% dtn) dt %>%
              data.table::setnames("Pos", "pos")

            if ("Core" %chin% dtn) dt %>% data.table::setnames("Core", "core")

          # fix netMHCpan allele output to match input
            if (command %like% "netMHCpan"){
              dt[, allele := allele %>%
                stringr::str_replace(fixed("*"), "")]
            }

          # set unique column names based on program
            data.table::setnames(dt, dt %>% names %exclude% "allele|nmer",
                                 paste0((dt %>% names %exclude% "allele|nmer"), "_", ptype))

          # name allele column for merge
            if (!"allele" %chin% (dt %>% names)) browser()
            dt %>%
              data.table::setnames("allele", ptype)

            return(dt)
                     })
   return(dtl)
 }



## ---- check_pred_tools
#' Internal function to check for netMHC tools and mhcflurry in `PATH`
#'
#' @export check_pred_tools
#' @md
check_pred_tools <- function(){

default_path <- paste0(system('echo $HOME', intern = TRUE),
                c(
                "/netMHC/netMHC-4.0",
                "/netMHC/netMHCII-2.2",
                "/netMHC/netMHCIIpan-3.1",
                "/netMHC/netMHCpan-3.0")) %>%
                paste(collapse = ":")

PATH_status <- list(
              mhcflurry = TRUE,
              netMHC = TRUE,
              netMHCpan = TRUE,
              netMHCII = TRUE,
              netMHCIIpan = TRUE)

Sys.setenv(PATH = paste0(default_path, ":", Sys.getenv("PATH")))

 if (suppressWarnings(system('which mhcflurry-predict 2> /dev/null', intern = TRUE)) %>%
        length == 0) {
        message("mhcflurry-predict is not in PATH\n       Download: https://github.com/hammerlab/mhcflurry")
      PATH_status$mhcflurry <- FALSE
      }
  if (suppressWarnings(system('which netMHC 2> /dev/null', intern = TRUE)) %>%
        length == 0) {
          message("netMHC is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHC/")
        PATH_status$netMHC <- FALSE
        }
  if (suppressWarnings(system('which netMHCpan 2> /dev/null', intern = TRUE)) %>%
        length == 0) {
          message("netMHCpan is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCpan/")
        PATH_status$netMHCpan <- FALSE
        }
  if (suppressWarnings(system('which netMHCII 2> /dev/null', intern = TRUE)) %>%
        length == 0) {
          message("netMHCII is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCII/")
        PATH_status$netMHCII <- FALSE
        }
  if (suppressWarnings(system('which netMHCIIpan 2> /dev/null', intern = TRUE)) %>%
        length == 0) {
          message("netMHCIIpan is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCIIpan/")
        PATH_status$netMHCIIpan <- FALSE
        }
        return(PATH_status)
}



## ---- run_netMHC
#' Internal function to run netMHC commands.
#'
#' @param dt Data table of prediction commands to run.
#'
#' @export run_netMHC
#' @md
run_netMHC <- function(dt){

  if (!"command" %chin% (dt %>% names))
    stop("dt must contain command column")

  # run commands
  esl <- parallel::mclapply(
         dt[, command],
          function(command){
          # run command
           es <- try(system(command, intern = TRUE))
          # if error, return empty dt
            if (es %>% class %>% .[1] == "try-error")
                return(data.table::data.table(status = 1, command = command))
            return(list(command, es))
            })

  dtl <- esl %>% collate_netMHC

      return(dtl)
}



## ---- write_nmers
#' Internal function to output nmers for MHC prediction to disk
#'
#' @param dt Data table of nmers.
#' @param type Character vector. Name of program to format for.
#'
#' @export write_nmers
#' @md

  write_nmers <- function(dt, type){
    if (dt %>% nrow == 0) return(NULL)
    if (!c("nmer", "nmer_l") %chin% (dt %>% names) %>% any)
          stop("dt must contain nmer and nmer_l columns")

    combs <- data.table::CJ(dt[, get(type)] %>% unique,
                            dt[, nmer_l] %>% unique)

    dto <- parallel::mclapply(1:nrow(combs), function(i) {

      dts <- dt[get(type) == combs$V1[i] & nmer_l == combs$V2[i]]

      # parallelize over 100 peptide chunks
      chunks <- ((dts %>% nrow)/100) %>% ceiling

      dto <- parallel::mclapply(dts %>% split(1:chunks), function(dtw){

        filename <- paste0(type, "_",
                    uuid::UUIDgenerate(), ".csv")


        # write out unique peptides for MHC type, length
        data.table::fwrite(dtw[, .(nmer)] %>% unique,
                          filename,
                          col.names = FALSE)

        return(data.table::data.table(
               type = type,
               allele = combs$V1[i],
               nmer_l = combs$V2[i],
               filename = filename))


        }) %>% data.table::rbindlist
      return(dto)

    }) %>% data.table::rbindlist

  return(dto)
  }



## ---- detect_hla
#' Replace HLA with matching type in netMHC format
#'
#' @param x Vector of HLA types named for program to convert to.
#' @param alleles Table of alleles to choose from.
#'
#' @export detect_hla
#' @md
 detect_hla <- function(x, alleles){

  prog <- deparse(substitute(x))

  for (hla in (x %>% unique)) {

    # match hla allele alelle (end|(not longer allele))
    hla_re <- paste0(hla, "($|[^0-9])")

    allele <- alleles[type == prog, allele %>%
      .[stringi::stri_detect_regex(., hla_re) %>% which]] %>%
      # match to first in case of multiple matches
      # e.g. netMHCII
      .[1]

      if (allele %>% length == 0) allele <- NA

      x[x == hla] <- allele
      }

    return(x)
  }



## ---- get_snpeff_annot
#' Internal function to extract SnpEff annotation information to a data table.
#'
#' @param dt Data table with INFO column.
#' @export get_snpeff_annot
#' @md
get_snpeff_annot <- function(dt){

    if (!"INFO" %chin% (dt %>% names)) stop("dt must contain INFO column")

    dt[, se := INFO %>%
      stringr::str_extract("ANN.*") %>%
      stringr::str_replace("ANN=[^\\|]+\\|", "")]

    # add a variant identifier
    suppressWarnings(dt[, snpeff_uuid :=
                  lapply(1:nrow(dt),
                  uuid::UUIDgenerate) %>% unlist])

    # abort if no variants passed filtering
    if (dt %>% nrow < 1) return(NULL)

    # spread SnpEff annotation over rows
    dt %>% tidyr::separate_rows("se", sep = ",")

    # extract info from snpeff annotation
      dt[, effect_type := se %>%
          stringr::str_extract("^[a-z0-9][^\\|]+")]
      dt[, ensembl_transcript_id := se %>%
          stringr::str_extract("(?<=\\|)(ENSMUST|ENST)[0-9]+")]
      dt[, ensembl_gene_id := se %>%
          stringr::str_extract("(?<=\\|)(ENSMUSG|ENSG)[0-9.]+(?=\\|)")]
      dt[, protein_change := se %>%
          stringr::str_extract("p\\.[^\\|]+")]
      dt[, cDNA_change := se %>%
          stringr::str_extract("c\\.[^\\|]+")]
      dt[, protein_coding := se %>%
          stringr::str_detect("protein_coding")]

      return(dt)

        }



## ---- get_metadata
#' Internal function to add metadata by `ensembl_transcript_id`
#'
#' @param dt Data table with `INFO` column.
#' @param humandb Character vector. One of `GRCh37` or `GRCh38`.
#' @param mousedb Character vector. One of `GRCm37` or `GRCm38`.
#'
#' @export get_metadata
#' @md
get_metadata <- function(dt,
                         humandb = "GRCh38",
                         mousedb = "GRCm38"){

  if (!"ensembl_transcript_id" %chin%
      (dt %>% names))
  stop("ensembl_transcript_id column missing")


  # set genome host
  if (!humandb %chin% c("GRCh37", "GRCh38")) stop("humandb set incorrectly")
  if (!mousedb %chin% c("GRCm37", "GRCm38")) stop("mousedb set incorrectly")
  if (humandb == "GRCh38") hhost <- "ensembl.org"
  if (humandb == "GRCh37") hhost <- "grch37.ensembl.org"
  if (mousedb == "GRCm38") mhost <- "ensembl.org"
  if (mousedb == "GRCm37") mhost <- "may2012.archive.ensembl.org"

    # remove version suffix
    dt[, ensembl_transcript_id :=
      ensembl_transcript_id %>%
      stringr::str_replace("\\.[0-9]$", "")]

    bmds <- vector()

    if (dt[, ensembl_transcript_id %>%
        stringr::str_detect("ENSMUST")] %>%
          stats::na.omit %>%
          unique %>%
          any) bmds <- c(bmds, "mmusculus_gene_ensembl")

    if (dt[, ensembl_transcript_id %>%
        stringr::str_detect("ENST")] %>%
          stats::na.omit %>%
          unique %>%
          any) bmds <- c(bmds, "hsapiens_gene_ensembl")

    message("Obtaining cDNA and peptide sequences using biomaRt")
    var_dt <- parallel::mclapply(bmds, function(i){

      if (i == "hsapiens_gene_ensembl") host <- hhost
      if (i == "mmusculus_gene_ensembl") host <- mhost

      mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                                 dataset = i,
                                                 host = host)

    if (i == "mmusculus_gene_ensembl") {
        trn <- dt[, ensembl_transcript_id %include% "ENSMUST" %>%
                  stats::na.omit %>%
                  unique]
              }

    if (i == "hsapiens_gene_ensembl") {
        trn <- dt[, ensembl_transcript_id %include% "ENST" %>%
                  stats::na.omit %>%
                  unique]
              }

   if (trn %>% length < 1) return(NULL)
   if (trn %>% length >= 1) {

    # obtain transcript metadata
      var_dt <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                    "external_gene_name", "ensembl_gene_id", "description", "chromosome_name",
                    "start_position", "end_position", "transcript_start", "transcript_end",
                    "transcript_length", "refseq_mrna"),
                         filters = c("ensembl_transcript_id"),
                         values = list(trn),
                         mart = mart) %>%
            data.table::as.data.table

      # obtain transcript cDNA and peptide sequences
      seqdtl <- parallel::mclapply(c("coding", "peptide"), function(j) {

         biomaRt::getSequence(type = "ensembl_transcript_id",
                     id = trn,
                     seqType = j,
                     mart = mart) %>% data.table::as.data.table
      })

      seqdt <- merge(seqdtl[[1]], seqdtl[[2]], by = "ensembl_transcript_id")
      var_dt <- merge(var_dt, seqdt, by = "ensembl_transcript_id")
      }

    return(var_dt)

    }) %>% data.table::rbindlist
    dt <- merge(dt, var_dt, by = "ensembl_transcript_id")

    return(dt)
}



## ---- make_cDNA
#' Internal function to create cDNA from HGVS notation
#'
#' @param dt Data table with INFO column.
#'
#' @export make_cDNA
#' @md
make_cDNA <- function(dt){

  if (!c("cDNA_type",
         "coding",
         "cDNA_locs",
         "cDNA_locl",
         "cDNA_seq"
         ) %chin%
      (dt %>% names) %>% all)
  stop("dt is missing columns")

  # initialize and set types

    dt[, coding_mut := coding]

    for (i in c("cDNA_locs", "cDNA_locl")) {
      set(dt, j = i, value = dt[, get(i) %>%
          as.integer])
    }

  # remove if nonsensical loc
    dt %<>% .[, coding_nchar := coding %>% nchar]

    dt %<>%
     .[!cDNA_locs > coding_nchar &
       !cDNA_locl > coding_nchar]

    dt[, coding_nchar := NULL]



  # paste substr around changed bases
  # substr is vectorized C, fast
  # mut sub nmers that are really wt
    # will be filtered after subpep generation


    # handle single base changes

      dt[cDNA_type == ">",
      coding_mut := coding_mut %>%
      # { controls . placement in paste
        { paste0(
           substr(., 1, cDNA_locs-1),
           cDNA_seq,
           substr(.,
                  cDNA_locs + 1,
                  nchar(.))) } ]

    # handle deletions
      # indices are inclusive
      # regexpr match for del

      dt[cDNA_type %like% "del",
      coding_mut := coding_mut %>%
      # { controls . placement in paste
        { paste0(
           substr(., 1, cDNA_locs-1),
           substr(.,
                  cDNA_locl + 1,
                  nchar(.))) } ]

    # handle insertions
      # regexpr match for del

      dt[cDNA_type %like% "ins",
      coding_mut := coding_mut %>%
      # { controls . placement in paste
        { paste0(
           substr(., 1, cDNA_locs),
           cDNA_seq,
           substr(.,
                  cDNA_locs + 1,
                  nchar(.))) } ]

      }



## ---- extract_cDNA
#' Internal function to extract cDNA changes from HGVS notation
#'
#' @param dt Data table with INFO column.
#'
#' @export extract_cDNA
#' @md
extract_cDNA <- function(dt){

  # check required cols

  if (!c("cDNA_change") %chin%
        (dt %>% names) %>% all)
    stop("cDNA_change missing")


  # extract HGVS DNA nomenclature
      # dups are just ins
      dt[, cDNA_change := cDNA_change %>%
            stringr::str_replace_all("dup", "ins")]

      dt[, cDNA_locs := cDNA_change %>%
            stringr::str_extract("[0-9]+") %>%
            as.integer]
      dt[, cDNA_locl := cDNA_change %>%
            stringr::str_extract("(?<=_)[0-9]+") %>%
            as.integer]
        # make cDNA_locl cDNA_locs for single bases
             dt[cDNA_locl %>% is.na, cDNA_locl := cDNA_locs]

      dt[, cDNA_type := cDNA_change %>%
            stringr::str_extract_all("[a-z]{3}|>") %>%
            unlist %>%
            paste(collapse = ""),
            by = 1:nrow(dt)]
      dt[, cDNA_seq := cDNA_change %>%
            stringr::str_extract("[A-Z]+$")]

  # filter out extraction errors
  # https://github.com/andrewrech/antigen.garnish/issues/3
    dt %<>% .[!cDNA_locs %>% is.na &
       !cDNA_locl %>% is.na &
       !cDNA_type %>% is.na]

  }



## ---- translate_cDNA
#' Internal function to translate cDNA to peptides
#'
#' @param v cDNA character vector without ambiguous bases.
#'
#' @export translate_cDNA
#' @md
  translate_cDNA <- function(v){

    parallel::mclapply(v, function(p){

     # protect vector length
       tryCatch({

          p %>%
            Biostrings::DNAString() %>%
            Biostrings::translate() %>%
                as.vector %>%
                as.character %>%
                paste(collapse = "")

        }, error = function(e) {
              return(NA)
              })

        }) %>% unlist
  }



## ---- garnish_variants
#' Intakes variants and returns an intersected data table for epitope prediction.
#'
#' Process variants annotated with [SnpEff](http://snpeff.sourceforge.net/). VCFs from matched samples are optionally intersected. [MuTect2](https://github.com/broadinstitute/gatk)/[Strelka2](https://github.com/Illumina/strelka)-derived VCFs are filtered for high confidence variants prior to intersection.
#'
#' @param vcfs Character vector. VFC files to import.
#'
#' @return A VCF as a data table with one unique SnpEff variant annotation per row, including:
#' * **sample_id**: sample identifier constructed from input \code{.bam} file names
#' * **se**: SnpEff annotation
#' * **effect_type**: SnpEff effect type
#' * **ensembl_transcript_id**: transcript effected
#' * **ensembl_gene_id**: gene effected
#' * **protein_change**: protein change in [HGVS](http://varnomen.hgvs.org/recommendations/DNA/) format
#' * **cDNA_change**: cDNA change in [HGVS](http://varnomen.hgvs.org/recommendations/protein/) format
#' * **protein_coding**: is the variant protein coding?
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download an example VCF
#'    dt <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'  # extract variants
#'    garnish_variants %T>%
#'    str
#'}
#' @export garnish_variants
#' @md
garnish_variants <- function(vcfs) {

  message("Loading VCFs")
  ivfdtl <- parallel::mclapply(vcfs %>% seq_along, function(ivf){

  # load dt
      vcf <-  vcfR::read.vcfR(vcfs[ivf], verbose = TRUE)

  # extract sample names from Mutect2 and Strelka command line for intersection

      sample_id <- vcf@meta %>%
                      stringr::str_extract_all("[^ ]+\\.bam") %>%
                      unlist %>%
                      unique %>%
                      basename %>%
                      sort %>%
                      paste(collapse = "_") %>%
                      stringr::str_replace("\\.bam", "")
      # extract vcf type
      vcf_type <- vcf@meta %>% unlist %>% stringr::str_extract(stringr::regex("Strelka|Mutect", ignore_case = TRUE)) %>% stats::na.omit %>%
                  unlist %>% data.table::first
      if (vcf_type %>% length == 0) vcf_type <- "other"

  # return a data table of variants

    vdt <- vcf@fix %>% data.table::as.data.table

    if (vcf@gt %>% length > 0) vdt <- cbind(vdt, vcf@gt %>% data.table::as.data.table)

    if(vdt %>% nrow < 1) return(data.table::data.table(sample_id = sample_id))

    # filter passing Strelka2 and muTect variants
    if(vcf_type == "Strelka") vdt <- vdt[FILTER == "PASS"]
    if(vcf_type == "Mutect") vdt <- vdt[INFO %>%
                                        stringr::str_extract("(?<=TLOD=)[0-9\\.]") %>%
                                        as.numeric > 6.0]
    vdt[, sample_id := sample_id]
    vdt[, vcf_type := vcf_type]

    # filter SnpEff warnings

    vdt %<>% get_snpeff_annot

    vdt %<>% .[!se %like% "ERROR_.*CHROMOSOME"]
    vdt %<>% .[!se %likef% "WARNING_SEQUENCE_NOT_AVAILABLE"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_INCOMPLETE"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS"]
    vdt %<>% .[!se %likef% "WARNING_TRANSCRIPT_NO_START_CODON"]

    return(vdt)
    })

    ivfdt <- ivfdtl %>% data.table::rbindlist

    merge_vcf <- function(dt, dt2){

    # a function to intersect annotated variants across VCFs using SnpEff

      sdt <- merge(dt, dt2[, .SD,
             .SDcols = c("CHROM",
              "POS", "REF", "cDNA_change")],
               by = c("CHROM", "POS", "REF",
              "cDNA_change")) %>% .[, vcf_type := "intersect"]
      return(sdt)

    }
  # return an intersected data table of variants
  sdt <- parallel::mclapply(ivfdt[, sample_id %>% unique], function(sn){

    # find data tables with matching sample names
    sdt <- lapply(ivfdtl, function(dt){

     dt[, sample_id %>% .[1]] == sn

    }) %>% unlist

  # merge all data tables with matching sample names
    if (ivfdtl[sdt] %>% length == 1) return(ivfdtl[[sdt %>% which]])
    if (ivfdtl[sdt] %>% length > 1) return(ivfdtl[sdt] %>% Reduce(merge_vcf, .))


  }) %>% data.table::rbindlist

  # select protein coding variants without NA
  sdt %<>% .[protein_coding == TRUE &
            !protein_change %>% is.na &
            !effect_type %>% is.na &
            effect_type %like% "insertion|deletion|missense|frameshift"]

  return(sdt)

}



## ---- garnish_predictions
#' Performs epitope prediction.
#'
#' Performs epitope prediction on a data table of missense mutations.
#'
#' @param dt Data table. Input data table from garnish_variants. Required columns: sample_id, ensembl_transcript_id (e.g. `ENST00000311936`), cDNA_change \href{http://varnomen.hgvs.org/recommendations/DNA/}{HGVS} notation (e.g. `c.718T>A`), space-separated MHC (e.g. `HLA-A*02:01 H-2-Kb HLA-DRB1*03:01`).
#' @param assemble Logical. Assemble data table?
#' @param generate Logical. Generate peptides?
#' @param predict Logical. Predict binding affinities?
#' @param humandb Character vector. One of "GRCh37" or "GRCh38".
#' @param mousedb Character vector. One of "GRCm37" or "GRCm38".
#' @return A data table of binding predictions including:
#' * **cDNA_seq**: mutant cDNA sequence
#' * **cDNA_locs**: starting index of mutant cDNA
#' * **cDNA_locl**: ending index of mutant cDNA
#' * **cDNA_type**: netMHC prediction tool output
#' * **frameshift**: frameshift variant?
#' * **coding**: wt cDNA sequence
#' * **coding_mut**: mutant cDNA sequence
#' * **pep_type**: type of peptide for this row
#' * **pep_mut**: mutant cDNA sequence
#' * **pep_wt**: wt cDNA sequence
#' * **mismatch_s**: starting index of mutant peptide sequence
#' * **mismatch_l**: ending index of mutant peptide sequence
#' * **mutant_loc**: index of mutant peptide for this row
#' * **nmer**: nmer for prediction
#' * **nmer_i**: index of nmer in sliding window
#' * ***_net**MHC MHCII ...]}: netMHC prediction tool output
#' * **mhcflurry_**[prediction ...]}: netMHC prediction tool output
#'
#' as well as a transcript description:
#' * description
#' * start_position
#' * end_position
#' * transcript_end
#' * transcript_length
#' * transcript_start
#' * peptide
#' * refseq_mrna
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download an example VCF
#'    dt <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'  # extract variants
#'    garnish_variants %>%
#'
#'  # predict neoepitopes
#'    garnish_predictions %T>%
#'    str
#'}
#'
#'\dontrun{
#'## ---- using an existing data table
#'
#'library(data.table)
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  dt <- data.table::data.table(
#'           sample_id = "test",
#'           ensembl_transcript_id =
#'           c("ENSMUST00000128119",
#'             "ENSMUST00000044250",
#'             "ENSMUST00000018743"),
#'           cDNA_change = c("c.4988C>T",
#'                           "c.1114T>G",
#'                           "c.718T>A"),
#'           MHC = c("HLA-A*02:01 HLA-DRB1*14:67",
#'                   "H-2-Kb H-2-IAd",
#'                   "HLA-A*01:47 HLA-DRB1*03:08")) %>%
#'  garnish_predictions %T>%
#'  str
#' }
#' @export garnish_predictions
#' @md
garnish_predictions <- function(dt,
                               assemble = TRUE,
                               generate = TRUE,
                               predict = TRUE,
                               humandb = "GRCh38",
                               mousedb = "GRCm38") {

  # remove temporary files on exit
  on.exit({
    message("Removing temporary files")
    list.files(pattern = "(netMHC|mhcflurry).*_[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12}\\.csv") %>% file.remove
  })

  dt %<>% data.table::as.data.table

  if (!c("sample_id", "ensembl_transcript_id", "cDNA_change", "MHC") %chin% (dt %>% names) %>% all) stop("Input dt must contain four columns:\n\nsample_id\nensembl_transcript_id (e.g. \"ENST00000311936\")\ncDNA_change in HGVS notation (e.g. \"c.718T>A\")\nMHC\n    e.g. \"HLA-A*02:01 HLA-A*03:01\" or\n         \"H-2-Kb H-2-Kb\" or \"HLA-DRB1*11:07 [second type]\"")

if (assemble){

    dt %<>% get_metadata(humandb = humandb, mousedb = mousedb)

    # extract cDNA changes and make cDNA
    dt %<>% extract_cDNA
    dt %<>% make_cDNA

    # translate protein sequences
    dt[, pep_mut := coding_mut %>% translate_cDNA]
    dt[, pep_wt := coding %>% translate_cDNA]

    dt[, cDNA_delta := ((coding_mut %>% nchar) - (coding %>% nchar)) / 3 ]

    dt[, frameshift := FALSE]

    # frameshifts have cDNA delta
    # not divisible by 3
      # does not handle stop codon loss
      # (this is acceptable because
      # no cDNA exists to know the readthrough)
    dt[cDNA_delta %% 3L != 0L, frameshift := TRUE]

    # remove variants with translated sequence-ensembl mismatch
    dt %<>% .[peptide == pep_wt]

    # remove stop codons

    for (i in dt %>% names %include% "^pep"){
      set(dt, j = i, value = dt[, get(i) %>%
                                stringr::str_extract_all("^[^\\*]+") %>%
                                unlist])

    }



## ---- create mutant peptide index

    # index first mismatch
    suppressWarnings({
    dt[, mismatch_s :=
        {
          (pep_wt %>%
               strsplit(split = "", fixed = TRUE) %>%
               unlist) !=
             (pep_mut %>%
               strsplit(split = "", fixed = TRUE) %>%
               unlist)
           } %>%
          which %>% .[1], by = 1:nrow(dt)]
          })

    # remove rows with matching transcripts
    dt %<>% .[!pep_wt == pep_mut]

    # initialize mismatch length
        dt[, mismatch_l := mismatch_s]

    # frameshifts are mutants until STOP
        dt[frameshift == TRUE,
        mismatch_l := pep_mut %>% nchar]

    # create mutant register for
    # non-frameshift insertions
        dt[mismatch_l > mismatch_s &
           frameshift == FALSE,
        mismatch_l := mismatch_s + (pep_mut %>% nchar) -
                      (pep_wt %>% nchar)]

    # deletions are mutants over site only
        dt[, mismatch_l := mismatch_s]

    # parallelized function to create a
    # space-separated character string
    # between two integers
      get_ss_str <- function(x, y) {
        mcMap(function(x, y) (x %>% as.integer):(y %>% as.integer) %>%
              paste(collapse = " "), x, y) %>% unlist
        }

    # create a space-separated vector of mutant indices
        dt[, mutant_loc := mismatch_s %>% as.character]
        dt[mismatch_l > mismatch_s, mutant_loc :=
            get_ss_str(mismatch_s, mismatch_l)]

    # warning if mutant_loc == NA
      if (dt[mutant_loc %>% is.na] %>%
            nrow > 1) {
        failn <- dt[mutant_loc %>% is.na] %>% nrow
        warning(paste0("Could not determine mutant index for ", failn, " records."))
      }
  }

if (generate) {
  message("Generating variants")
  # generation a uuid for each unique variant
  suppressWarnings(dt[, var_uuid :=
                  parallel::mclapply(1:nrow(dt),
                  uuid::UUIDgenerate) %>% unlist])

  # separate over mutant indices

    if (dt[, mutant_loc %>% unique] %>%
        stringr::str_detect(" ") %>% any) {
      dts <- dt %>% tidyr::separate_rows("mutant_loc", sep = " ")
    } else {
      dts <- dt
    }

  # convert back to numeric
  dts[, mutant_loc := mutant_loc %>% as.numeric]

  # generate a data table of unique variants for peptide generation

    dtnfs <- dts[frameshift == FALSE]
    dtfs <- dts[frameshift == TRUE]

    basepep_dt <- data.table::rbindlist(list(
         # take pep_wt for non-fs for DAI calculation
           data.table::rbindlist(list(
                dtnfs %>%
                data.table::copy %>%
                .[, pep_base := pep_wt] %>%
                .[, pep_type := "wt"],
                dtnfs %>%
                data.table::copy %>%
                .[, pep_base := pep_mut] %>%
                .[, pep_type := "mutnfs"]
                 )),
         # take only pep_mut for fs
            dtfs %>%
            data.table::copy %>%
            .[, pep_base := pep_mut] %>%
            .[, pep_type := "mutfs"])) %>%

         # take unique peptides
            .[, .SD, .SDcols = c("var_uuid",
                                 "pep_type",
                                 "pep_base",
                                 "mutant_loc")] %>%
            unique

    if (basepep_dt %>% nrow == 0)
      return("no variants for peptide generation")

    sink(file = "/dev/null")
    nmer_dt <- garnish_predictions_worker(basepep_dt) %>% .[, nmer_l := nmer %>% nchar]
    sink()

     dt <- merge(dt, nmer_dt,
        by = "var_uuid",
        all.x = TRUE)

    # generation a uuid for each unique nmer
    suppressWarnings(dt[, nmer_uuid :=
                    parallel::mclapply(1:nrow(dt),
                    uuid::UUIDgenerate) %>% unlist])

     dt %<>% get_DAI_uuid

    # remove mut == wt by sample_id
    # remove wt missing mut counterpart

    for (id in dt[, sample_id %>% unique]) {

     # get wt nmers for fast !chmatch
     id_wt_nmers <- dt[sample_id == id &
                        pep_type == "wt",
                        nmer %>% unique]

     dt %<>% .[pep_type == "wt" |
               sample_id != id |
               (sample_id == id &
                pep_type != "wt" &
                !nmer %chin% id_wt_nmers)]
         }

    # remove wt without a matched mut
     unmatched_dai <- dt[, .N, by = dai_uuid] %>% . [N == 1, dai_uuid]
     if (unmatched_dai %>% length > 0) dt[!dai_uuid %chin% unmatched_dai]

}

if (predict) {

    dtl <- dt %>% get_pred_commands

  # run commands
   if (!check_pred_tools() %>% unlist %>% all) {
    warning("Missing prediction tools in PATH, returning data table without predictions.")
    return(dtl[[1]])
  }

    # netMHC
    message("Running netMHC in parallel")
    dto <- run_netMHC(dtl[[2]])

    # mhcflurry
    message("Running mhcflurry in parallel")
    list.files(pattern = "mhcflurry_input.*csv") %>%

    mclapply(., function(x){
      paste0("mhcflurry-predict ", x, " > ", x %>%
            stringr::str_replace("input", "output")) %>%
      system
          })

   dt <- merge_predictions(dto, dtl[[1]])

}
   return(dt)
}



## ---- garnish_predictions_worker
#' Internal function for parallelized `nmer` creation.
#'
#' @param dt Data table. Input data table from `garnish_predictions`.
#'
#' @export garnish_predictions_worker
#' @md
garnish_predictions_worker <- function(dt) {


if (!c("var_uuid",
       "pep_type",
       "pep_base",
       "mutant_loc") %chin% (dt %>% names) %>% all)
  stop("dt is missing columns")

  dt %<>% data.table::as.data.table

  lines <- nrow(dt)

  # a function to generate sliding
  # window n-mers over an index of interest

  message("Generating nmers")
  nmer_dt <- parallel::mclapply(1:nrow(dt),
                                function(n){


      ## --- Write peptide fragments

        # for every peptide length
        nmer_dt <- lapply((15:8), function(pl) {

          mut_frag_t <- dt$pep_base[n] %>% strsplit("",
                              fixed = TRUE) %>% unlist
          mut_frag_loc <- dt$mutant_loc[n]

          # if the peptide is not long enough, return
          if (!(mut_frag_t %>% length) >= pl) return(NULL)

            # re-register peptide if the mutant index
            # is not centered due to back truncation

            # trim peptide front to mutant site
            if (mut_frag_loc > pl){
              rml <- mut_frag_loc - pl + 1
              mut_frag_t <- mut_frag_t[rml:(length(mut_frag_t))]
              mut_frag_loc <- pl
              }

            # trim peptide back to mutant site
            fpl <- mut_frag_loc + pl - 1
            if (fpl < length(mut_frag_t))
                mut_frag_t <- mut_frag_t[1:(fpl)]

          # slide across the peptide window
          # create (n = pl )-mers wrapped in
          # sync to prevent zoo::rollapply stdout

          nmers <- zoo::rollapply(mut_frag_t, pl,
                        print,
                        partial = FALSE,
                        align = "left") %>%
          apply(1, function(pmr){
          paste(pmr, sep = "", collapse = "")
          })


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


    }) %>% data.table::rbindlist %>% unique

      return(nmer_dt)


    }) %>% data.table::rbindlist %>% unique

  return(nmer_dt)
}



## ---- garnish_summary
#' Summarize epitope prediction.
#'
#' Calculate neoepitope summary statistics over samples.
#'
#' @param dt Data table. Prediction output from `garnish_predictions`.
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download an example VCF
#'    dt <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'  # extract variants
#'    garnish_variants %>%
#'
#'  # predict neoepitopes
#'    garnish_predictions %>%
#'
#'  # summarize predictions
#'    garnish_summary %T>%
#'    print
#'}
#' @return A summary data table of \code{dt} by \code{sample_id} with the following columns:
#' * **priority_neos**: neoepitopes that meet both classic and alternative criteria, or that meet classic criteria and are derived from frameshift mutations
#' * **classic_neos**: high affinity neoepitopes
#' * **classic_top_score**: sum of the top three affinity scores (1 / MHC affinity (nM))
#' * **alt_neos**: mutant nmers predicted to bind MHC with greatly improved affinity relative to their non-mutated counterpart
#' * **alt_top_score**: sum of the top three `nmer` differential agretopicity indices; differential agretopicity index (DAI) is the ratio of MHC binding affinity between mutant and wt peptide.
#' * **nmers**: mutant nmers
#' * **variants**: genomic variants evaluated
#' * **transcripts**: transcripts evaluated
#' * **predictions**: wt and mutant predictions performed
#' * **mhc_binders**: nmers predicted to at least minimally bind MHC
#' @export garnish_summary
#' @md
garnish_summary <- function(dt){

# summarize over unique nmers

  dt %<>% data.table::as.data.table %>%
    data.table::copy %>%
    unique(by = c("pep_type",
                  "MHC",
                  "nmer"))

  dt <- dt[DAI != Inf & DAI != -Inf & Consensus_scores != Inf & Consensus_scores != -Inf]

    # function to sum top values of a numeric vector
      sum_top_v <- function(x, value = 3){

        x %<>% sort %>% rev
        return(sum(x[1:value]))
      }

  dt %>% data.table::setkey(sample_id)

  dtn <- parallel::mclapply(dt[, sample_id %>% unique], function(id){

    dt <- dt[sample_id == id]

      return(
          data.table::data.table(
          sample_id = id,
          priority_neos = dt[DAI > 10 & (
                       Consensus_scores < 50 |
                        (frameshift == TRUE &
                         (
                          mutant_loc > (mismatch_s + 2)
                          )
                          )
                       )] %>% nrow,
          classic_neos = dt[Consensus_scores < 50] %>% nrow,
          classic_top_score = dt[Consensus_scores < 5000, (1/Consensus_scores) %>% sum_top_v],
          alt_neos = dt[Consensus_scores < 5000 & DAI > 10] %>% nrow,
          alt_top_score = dt[Consensus_scores < 5000, DAI %>% sum_top_v],
          mhc_binders = dt[Consensus_scores < 5000] %>% nrow,
          variants = dt[, var_uuid %>% unique] %>% length,
          transcripts = dt[ensembl_transcript_id %>% unique] %>% nrow,
          predictions = dt[pep_type %like% "mut"] %>% nrow,
          nmers = dt[pep_type %like% "mut", nmer %>% unique] %>% length
          ))

    }) %>% data.table::rbindlist

  return(dtn)
}
