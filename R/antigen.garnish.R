


## ---- get_DAI_uuid
#' Pair peptides from missense sites by a common UUID for DAi calculations
#'
#' @param dt Data table of nmers.
#'
#' @export get_DAI_uuid

get_DAI_uuid <- function(dt){

  ### package building
  .sinew <- function(){


      invisible({
        testthat::compare()
        dt.inflix::allduplicated()
        colorspace::RGB()
                })
       ###
  }


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


## ---- check_pred_tools
#' Check for netMHC tools and mhcflurry in PATH
#'
#' @export check_pred_tools
#'

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
#' Run netMHC commands.
#'
#' @param dt Data table of commands to run.
#'
#' @export run_netMHC

run_netMHC <- function(dt){

  if (!"command" %chin% (dt %>% names))
    stop("dt must contain command column")

   if (!check_pred_tools() %>% unlist %>% all) {
    stop("Missing prediction tools in PATH")
  }

  dtl <- parallel::mclapply(
         dt[, command],
          function(command){
          # run command
           es <- try(system(command, intern = TRUE))
          # if error, return empty dt
            if (es %>% class %>% .[1] == "try-error") {
              return(data.table::data.table(status = 1, command = command))
            }

          # parse results
            # isolate table and header

              dtl <- es %exclude%
                "^\\#|----|Number|Distance|threshold|version|^$" %>%
                stringr::str_replace("^[ ]+", "")
              dtn <- dtl %include%
                "[Pp]eptide" %>%
                strsplit("[ ]+") %>%
                unlist %exclude% "^$|^Bind|^Level"

                if (!dtn %>% is.character) browser()

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

          # make netMHC names consistent
            if ("Peptide" %chin% (dt %>% names)) dt %>% data.table::setnames("Peptide", "peptide")
            if ("peptide" %chin% (dt %>% names)) dt %>% data.table::setnames("peptide", "nmer")
            if ("Affinity(nM)" %chin% (dt %>% names)) dt %>% data.table::setnames("Affinity(nM)", "affinity(nM)")
            if ("Aff(nM)" %chin% (dt %>% names)) dt %>% data.table::setnames("Aff(nM)", "affinity(nM)")
            if ("HLA" %chin% (dt %>% names)) dt %>%
              data.table::setnames("HLA", "allele")
            if ("Allele" %chin% (dt %>% names)) dt %>% data.table::setnames("Allele", "allele")
            if ("Pos" %chin% (dt %>% names)) dt %>%
              data.table::setnames("Pos", "pos")
            if ("Icore" %chin% (dt %>% names)) dt %>% data.table::setnames("Icore", "icore")
            if ("iCore" %chin% (dt %>% names)) dt %>% data.table::setnames("iCore", "icore")
            if ("Core" %chin% (dt %>% names)) dt %>% data.table::setnames("Core", "core")

          # fix netMHCpan allele output to match input
            if (command %like% "netMHCpan"){
              dt[, allele := allele %>%
                stringr::str_replace(fixed("*"), "")]
            }

          # set unique column names based on program
            data.table::setnames(dt, dt %>% names %exclude% "allele|nmer",
                                 paste0((dt %>% names %exclude% "allele|nmer"), "_", ptype))

          # name allele column for merge
            dt %>%
              data.table::setnames("allele", ptype)

            return(dt)
                     })
      return(dtl)
}



## ---- write_nmers
#' Write MHC prediction nmers to disk
#'
#' @param dt Data table of nmers.
#' @param type Character vector. Name of program to format for.
#'
#' @export write_nmers


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

 detect_hla <- function(x, alleles){

  prog <- deparse(substitute(x))

  for (hla in (x %>% unique)) {

    allele <- alleles[type == prog, allele %>%
      .[stringi::stri_detect_fixed(., hla) %>% which]]

      if (allele %>% length == 0) allele <- NA

      x %<>% stringr::str_replace_all(stringr::fixed(hla), allele)
      }

    return(x)
  }



## ---- get_snpeff_annot
#' Extract snpeff annotation information to a data table.
#'
#' @param dt Data table with INFO column.
#' @export get_snpeff_annot

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



## ---- aa_convert
#' Convert amino acid codes quickly
#'
#' @param x Character vector to convert.
#'
#' @export aa_convert

aa_convert <- function(x){

if (x %>% stats::na.omit %>%
    is.null) return(x)
if (x %>% stats::na.omit %>%
    length == 0) return(x)

if (x %>% stats::na.omit %>%
    stringr::str_detect("[A-Z][a-z]{2}") %>% any){
     x %<>% stringr::str_replace(stringr::fixed("Ala"), "A") %>%
            stringr::str_replace(stringr::fixed("Arg"), "R") %>%
            stringr::str_replace(stringr::fixed("Asn"), "N") %>%
            stringr::str_replace(stringr::fixed("Asp"), "D") %>%
            stringr::str_replace(stringr::fixed("Cys"), "C") %>%
            stringr::str_replace(stringr::fixed("Glu"), "E") %>%
            stringr::str_replace(stringr::fixed("Gln"), "Q") %>%
            stringr::str_replace(stringr::fixed("Gly"), "G") %>%
            stringr::str_replace(stringr::fixed("His"), "H") %>%
            stringr::str_replace(stringr::fixed("Ile"), "I") %>%
            stringr::str_replace(stringr::fixed("Leu"), "L") %>%
            stringr::str_replace(stringr::fixed("Lys"), "K") %>%
            stringr::str_replace(stringr::fixed("Met"), "M") %>%
            stringr::str_replace(stringr::fixed("Phe"), "F") %>%
            stringr::str_replace(stringr::fixed("Pro"), "P") %>%
            stringr::str_replace(stringr::fixed("Ser"), "S") %>%
            stringr::str_replace(stringr::fixed("Thr"), "T") %>%
            stringr::str_replace(stringr::fixed("Trp"), "W") %>%
            stringr::str_replace(stringr::fixed("Tyr"), "Y") %>%
            stringr::str_replace(stringr::fixed("Val"), "V")

} else {
     x %<>% stringr::str_replace(stringr::fixed("A"), "Ala") %>%
            stringr::str_replace(stringr::fixed("R"), "Arg") %>%
            stringr::str_replace(stringr::fixed("N"), "Asn") %>%
            stringr::str_replace(stringr::fixed("D"), "Asp") %>%
            stringr::str_replace(stringr::fixed("C"), "Cys") %>%
            stringr::str_replace(stringr::fixed("E"), "Glu") %>%
            stringr::str_replace(stringr::fixed("Q"), "Gln") %>%
            stringr::str_replace(stringr::fixed("G"), "Gly") %>%
            stringr::str_replace(stringr::fixed("H"), "His") %>%
            stringr::str_replace(stringr::fixed("I"), "Ile") %>%
            stringr::str_replace(stringr::fixed("L"), "Leu") %>%
            stringr::str_replace(stringr::fixed("K"), "Lys") %>%
            stringr::str_replace(stringr::fixed("M"), "Met") %>%
            stringr::str_replace(stringr::fixed("F"), "Phe") %>%
            stringr::str_replace(stringr::fixed("P"), "Pro") %>%
            stringr::str_replace(stringr::fixed("S"), "Ser") %>%
            stringr::str_replace(stringr::fixed("T"), "Thr") %>%
            stringr::str_replace(stringr::fixed("W"), "Trp") %>%
            stringr::str_replace(stringr::fixed("Y"), "Tyr") %>%
            stringr::str_replace(stringr::fixed("V"), "Val")
}
return(x)

}



## ---- get_metadata
#' Add metadata using ensembl_transcript_ID
#'
#' @param dt Data table with INFO column.
#'
#' @export get_metadata

get_metadata <- function(dt){

  if (!"ensembl_transcript_id" %chin%
      (dt %>% names))
  stop("ensembl_transcript_id column missing")

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
    var_dt <- lapply(bmds, function(i){

      mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                                 dataset = i,
                                                 host = 'ensembl.org')

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
#' Create cDNA from hgvs nomenclature (http://varnomen.hgvs.org/)
#'
#' @param dt Data table with INFO column.
#'
#' @export make_cDNA

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
#' Extract cDNA changes from hgvs nomenclature (http://varnomen.hgvs.org/)
#'
#' @param dt Data table with INFO column.
#'
#' @export extract_cDNA

extract_cDNA <- function(dt){

  # check required cols

  if (!c("cDNA_change") %chin%
        (dt %>% names) %>% all)
    stop("cDNA_change missing")


  # extract hgvs DNA nomenclature
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

## ---- ftrans
#' Translate cDNA to peptide fast
#'
#' @param v cDNA character vector without ambiguous bases.
#'
#' @export ftrans

  ftrans <- function(v){

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
#' Process raw variants from a \href{https://github.com/broadinstitute/gatk}{MuTect2}/\href{https://github.com/Illumina/strelka}{Strelka2} - \href{https://github.com/pcingola/SnpEff}{SnpEff} variant annotation pipeline and filters for neoepitope prediction. Hg38 (human) and GRCm38 (murine) variant calls are required. Mutect2 and Strelka variant threshold prior to intersection were empirically established to limit false positives.
#'
#' @param vcfs Character vector. VFC files to import.
#'
#' @return A data table of intersected protein coding variants for neoepitope prediction.
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
#'    garnish_variants %T>%
#'    str
#'
#'  # example output
#'    dt <- data.table::fread(
#'     "http://get.rech.io/antigen.garnish_example_output.txt") %T>%
#'     str
#'}
#' @export garnish_variants

garnish_variants <- function(vcfs) {

  message("Loading VCFs")
  ivfdtl <- lapply(vcfs %>% seq_along, function(ivf){

  # load dt
      vcf <-  vcfR::read.vcfR(vcfs[ivf], verbose = TRUE)

  # extract sample names from Mutect2 and Strelka command line for intersection

      sample_id <- vcf@meta %>% stringr::str_extract_all("[^ ]+\\.bam") %>% unlist %>% basename %>% paste(collapse = "_") %>% stringr::str_replace("\\.bam", "")
      # extract vcf type
      vcf_type <- vcf@meta %>% unlist %>% stringr::str_extract(stringr::regex("Sdtrelka|Mudtect", ignore_case = TRUE)) %>% stats::na.omit %>%
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

    vdt %<>% get_snpeff_annot

    return(vdt)
    })

    ivfdt <- ivfdtl %>% data.table::rbindlist

    merge_vcf <- function(dt, dt2){

    # a function to intersect annotated variants across VCFs using SnpEff

      sdt <- merge(dt, dt2[, .SD,
             .SDcols = c("CHROM",
              "POS",
              "REF",
              "cDNA_change")],
               by = c("CHROM",
              "POS",
              "REF",
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
#' @param dt Data table. Input data table from garnish_variants. \href{http://get.rech.io/antigen.garnish_example_input.txt}{Example.}
#' @param assemble Logical. Assemble data table?
#' @param generate Logical. Generate peptides?
#' @param predict Logical. Predict binding affinities?
#'
#' @return A data table of neoepitopes.
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
#'    garnish_predictions %T>%
#'    str
#'}
#'
#' @export garnish_predictions

garnish_predictions <- function(dt,
                               assemble = TRUE,
                               generate = TRUE,
                               predict = TRUE) {

  if (!"data.frame" %chin% (dt %>% class)) stop("dt must be a data frame or data table")

  if (!c("sample_id", "ensembl_transcript_id", "cDNA_change", "MHC") %chin% (dt %>% names) %>% all) stop("Input dt must contain four columns:\n\nsample_id\nensembl_transcript_id\ncDNA_change HGVS notation)\nMHC\n    e.g. 'HLA-A*02:01 HLA-A*03:01' or\n    'H-2-Kb H-2-Kb' or HLA-DRB1*11:07 [second type]'\n\n See: http://get.rech.io/antigen.garnish_example_input.txt")

  dt %<>% data.table::as.data.table

if (assemble){

    dt %<>% get_metadata

    dt %<>% extract_cDNA

    dt[, frameshift := FALSE]
    dt[protein_change %like% "fs$", frameshift := TRUE]

    dt %<>% make_cDNA

    dt[, pep_mut := coding_mut %>% ftrans]
    dt[, pep_wt := coding %>% ftrans]

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
    dt[, mismatch_f :=
        {
          (pep_wt %>%
               strsplit(split = "", fixed = TRUE) %>%
               unlist) !=
             (pep_mut %>%
               strsplit(split = "", fixed = TRUE) %>%
               unlist)
           } %>%
          which %>% .[1], by = 1:nrow(dt)]

    # initialize mismatch length
        dt[, mismatch_l := mismatch_f]

    # frameshifts are mutants until STOP
        dt[frameshift == TRUE,
        mismatch_l := pep_mut %>% nchar]

    # create mutant register for
    # non-frameshift insertions
        dt[mismatch_l > mismatch_f &
           frameshift == FALSE,
        mismatch_l := mismatch_f + (pep_mut %>% nchar) -
                      (pep_wt %>% nchar)]

    # deletions are mutants over site only
        dt[, mismatch_l := mismatch_f]

    # parallelized function to create a
    # space-separated character string
    # between two integers
      get_ss_str <- function(x, y) {
        mcMap(function(x, y) (x %>% as.integer):(y %>% as.integer) %>%
              paste(collapse = " "), x, y) %>%
              unlist
        }

    # create a space-separated vector of mutant locations
        dt[, mutant_loc := mismatch_f %>% as.character]
        dt[mismatch_l > mismatch_f, mutant_loc :=
            get_ss_str(mismatch_f, mismatch_l)]

  }

if (generate) {
  message("Generating variants")
  # generation a uuid for each unique variant
  suppressWarnings(dt[, var_uuid :=
                  parallel::mclapply(1:nrow(dt),
                  uuid::UUIDgenerate) %>% unlist])

  # separate over mutant locations

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

  if (dt[, MHC %>% unique] %>%
      stringr::str_detect(" ") %>% any) dt %<>%
      tidyr::separate_rows("MHC", sep = " ")

  # get available MHC alleles for predictions

  alleles <- data.table::rbindlist(
                         list(
  system.file("extdata",
      "netMHC_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHC"],
  system.file("extdata",
      "netMHCpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCpan"],
  system.file("extdata",
      "mhcflurry_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcflurry"],
  system.file("extdata",
      "netMHCII_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCII"],
  system.file("extdata",
      "netMHCIIpan_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE) %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "netMHCIIpan"]
                      ))

  # generate csv for mhcflurry predictions

    dt[MHC %chin% alleles[type == "mhcflurry", allele] &
    nmer_l < 15,
          .SD, .SDcols = c("MHC", "nmer")] %>%
    data.table::copy %>%
      data.table::setnames(c("MHC", "nmer"), c("allele", "peptide")) %>%
      unique %>%
      data.table::fwrite("mhcflurry_input.csv")

  # generate matchable MHC substring for netMHC tools

    dt[, netMHCpan := MHC %>% stringr::str_replace_all(stringr::fixed("*"), "")]
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

  # run commands

    # netMHC
    message("Running netMHC")
    ag_out_raw <- run_netMHC(dtfn)
    saveRDS(ag_out_raw, "ag_out_raw.RDS")

    # mhcflurry
    message("Running mhcflurry")
    system("mhcflurry-predict mhcflurry_input.csv > mhcflurry_output.csv")

    message("Merging output")
  # merge netMHC by program type
    progl <- lapply(ag_out_raw %>% seq_along, function(dti) {
     ag_out_raw[[dti]]$command[1] %>% stringr::str_extract("net[A-Za-z]+")
     })

  # merge netMHC output
    netmprogs <- progl %>% unique %>% unlist

    for (ptype in netmprogs) {

      dt <- merge(dt, ag_out_raw[(progl == ptype) %>% which] %>%
            data.table::rbindlist, by = c("nmer", ptype), all.x = TRUE)
    }

  # read and merge mhcflurry output
   fdt <- data.table::fread("mhcflurry_output.csv") %>%
      data.table::setnames(c("allele", "peptide"), c("MHC", "nmer"))

   dt <- merge(dt, fdt, by = c("nmer", "MHC"), all.x = TRUE)

  # calculate netMHC consensus score, preferring non-*net tools
      for (col in (dt %>% names %include% "aff|[Rr]ank|Consensus_scores")) {
        suppressWarnings(set(dt, j = col, value = dt[, get(col) %>% as.numeric]))
      }

    # get vector of netMHC scores
      cols <- dt %>% names %includef% c("affinity(nM)")

      dt[, Consensus_scores := c(.(get(cols)) %>%
      stats::na.omit, NA) %>% data.table::first, by = 1:nrow(dt)]

  # take average of mhcflurry and best available netMHC tool
    dt[, Consensus_scores := c(Consensus_scores,
                                   mhcflurry_prediction) %>%
                                  mean(na.omit = TRUE),
                                  by = 1:nrow(dt)]
  # remove nmers without predictions

    dt <- dt[!Consensus_scores %>% is.na]

  # calculate DAI

    dt[, DAI := NA %>% as.numeric]

    data.table::setkey(dt, pep_type, dai_uuid)

    # dai_uuid is always length 2
    # calculate DAI by dividing

    dt[!dai_uuid %>% is.na,
    DAI := Consensus_scores[2] /
      Consensus_scores[1], by = dai_uuid]

    return(dt)
}
}


## ---- garnish_predictions_worker
#' Parallelized worker function for garnish_predictions
#'
#' @param dt Data table. Input data table from garnish_predictions.
#'
#' @return A data table summary of neoepitope by sample_id.
#'
#' @export garnish_predictions_worker

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

    tryCatch({

      ## --- Write peptide fragments

        # for every peptide length
        nmer_dt <- lapply((15:8), function(pl) {

          mut_frag_t <- dt$pep_base[n] %>% strsplit("",
                              fixed = TRUE) %>% unlist
          mut_frag_loc <- dt$mutant_loc[n]

          # if the peptide is not long enough, return
          if (!(mut_frag_t %>% length) >= pl) return(NULL)

            # re-register peptide if the mutant location
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

    }, error = function(e) {
          print(cat("ERROR :", conditionMessage(e), "\n"))
          return(NULL)
          })

    }) %>% data.table::rbindlist %>% unique

  return(nmer_dt)
}



## ---- garnish_summary
#' Summarize epitope prediction.
#'
#' Calculate neoepitope summary statistics over samples.
#'
#' @param dt Data table. Prediction output from garnish_predictions.
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
#'
#'  # example output
#'    dt <- data.table::fread(
#'     "http://get.rech.io/antigen.garnish_example_summary.txt") %T>%
#'     str
#'}
#'
#' @export garnish_summary

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
          variants = dt[, var_uuid %>% unique] %>% length,
          priority_neos = dt[Consensus_scores < 50 & DAI > 10] %>% nrow,
          classic_neos = dt[Consensus_scores < 50] %>% nrow,
          alt_neos = dt[Consensus_scores < 5000 & DAI > 10] %>% nrow,
          alt_neos_top = dt[Consensus_scores < 5000, DAI %>% sum_top_v],
          classic_neos_top = dt[Consensus_scores < 5000, (1/Consensus_scores) %>% sum_top_v],
          binders = dt[Consensus_scores < 5000] %>% nrow,
          nmers = dt[pep_type %like% "mut", nmer %>% unique] %>% length,
          predictions = dt[pep_type %like% "mut"] %>% nrow
          ))

    }) %>% data.table::rbindlist

  return(dtn)
}
