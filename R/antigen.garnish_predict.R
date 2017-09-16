


## ---- make_DAI_uuid
#' Internal function to pair peptides from missense sites by a common UUID for DAI calculations
#'
#' @param dt Data table of nmers.
#'
#' @export make_DAI_uuid
#' @md

make_DAI_uuid <- function(dt){

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
#' Internal function to merge input data table and prediction results
#'
#' @param l Output list from run_netMHC
#' @param dt Input data table.
#' @export merge_predictions
#' @md

merge_predictions <- function(l, dt){

      message("Merging output")

      # merge netMHC by program type
        progl <- lapply(l %>% seq_along, function(dti){
          l[[dti]]$command[1] %>% stringr::str_extract("net[A-Za-z]+")
          })

       # merge netMHC output

        for (ptype in (progl %>% unique %>% unlist)){

          dt <- merge(dt, l[(progl == ptype) %>% which] %>%
                data.table::rbindlist, by = c("nmer", ptype), all.x = TRUE,
                allow.cartesian = TRUE)
        }
          dt %<>% unique

      # read and merge mhcflurry output

        f_flurry <- list.files(pattern = "mhcflurry_output.*csv")

        if (f_flurry %>% length > 0){
              fdt <- lapply(f_flurry, function(x){
                  if (data.table::fread(x) %>% nrow > 0){
                    return(data.table::fread(x))
                  } else {
                    return(NULL)
                  }
                }) %>% data.table::rbindlist %>%
                    data.table::setnames(c("allele", "peptide"), c("MHC", "nmer"))
                    dt <- merge(dt, fdt %>% unique, by = c("nmer", "MHC"), all.x = TRUE)
                    dt %<>% unique
                  }

      # read and merge mhcnuggets output

        f_nuggets <- list.files(pattern = "mhcnuggets_output.*csv")

        if (f_nuggets %>% length > 0){
                  nugdt <- lapply(f_nuggets, function(x){
                          if (data.table::fread(x) %>% nrow >0 ){
                            return(
                             data.table::fread(x) %>%
                              .[, mhcnuggets := basename(x) %>%
                                stringr::str_extract(pattern = "(?<=_)(H-2-.*(?=_))|(HLA).*(?=_)")] %>%
                              .[, tool := basename(x) %>%
                                stringr::str_extract(pattern = "(gru)|(lstm)")])
                           } else {
                              return(NULL)
                              }
                        }) %>%

                        data.table::rbindlist(fill = TRUE) %>%
                        data.table::setnames(c("Building", "model"),
                                             c("nmer", "mhcnuggets_prediction")) %>%
                        .[tool == "gru", mhcnuggets_pred_gru := mhcnuggets_prediction] %>%
                        .[tool == "lstm", mhcnuggets_pred_lstm := mhcnuggets_prediction] %>%
                        .[, c("nmer", "mhcnuggets", "mhcnuggets_pred_gru", "mhcnuggets_pred_lstm")]

                 nugdt <- merge(
                             nugdt[!is.na(mhcnuggets_pred_lstm),
                                    .SD,
                                    .SDcols = c("nmer",
                                    "mhcnuggets",
                                    "mhcnuggets_pred_lstm")] %>% unique,
                             nugdt[!is.na(mhcnuggets_pred_gru),
                                    .SD,
                                    .SDcols = c("nmer",
                                    "mhcnuggets",
                                    "mhcnuggets_pred_gru")] %>% unique,
                                by = c("nmer", "mhcnuggets"),
                                all = TRUE)

                dt <- merge(dt, nugdt %>%
                            unique,
                            by = c("nmer", "mhcnuggets"),
                            all.x = TRUE)
              }
        dt %<>% unique

      # calculate netMHC consensus score, preferring non-*net tools
         for (col in (dt %>% names %include% "aff|[Rr]ank|Consensus_scores")){
              suppressWarnings({
              set(dt, j = col, value = dt[, get(col) %>% as.numeric])
              })
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
             .[, .(best_netMHC =
        # define Consensus_score
                value %>%
                na.omit %>%
                .[1]), by = c("nmer", "MHC")] %>%
             .[!best_netMHC %>% is.na]

        # merge back
        dt %<>% merge(dtm, by = c("nmer", "MHC"))

        # take average of mhcflurry, mhcnuggets, and best available netMHC tool
        cols <- cols <- dt %>% names %include% "(best_netMHC)|(mhcflurry_prediction$)|(mhcnuggets_pred_gru)|(mhcnuggets_pred_lstm)"
        dt[, Consensus_scores :=  mean(as.numeric(.SD), na.rm = TRUE),
                                        by = 1:nrow(dt), .SDcols = cols]

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

  # replace all with all MHC types
    dt[MHC == "all", MHC :=
    system.file("extdata",
      "all_alleles.txt", package = "antigen.garnish") %>%
       data.table::fread(header = FALSE, sep = "\t") %>%
       .$V1 %>%
       paste(collapse = " ")
  ]
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
                      .[, type := "netMHCIIpan"],
  system.file("extdata",
       "mhcnuggets_gru_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcnuggets_gru"],
  system.file("extdata",
      "mhcnuggets_lstm_alleles.txt", package = "antigen.garnish") %>%
                      data.table::fread(header = FALSE, sep = "\t") %>%
                      data.table::setnames("V1", "allele") %>%
                      .[, type := "mhcnuggets_lstm"]
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
            "mhcflurry_input_",  uuid::UUIDgenerate() %>% substr(1, 18), ".csv"))
    }

  # generate input for mhcnuggets predictions

    dt[, mhcnuggets := toupper(MHC) %>%
       stringr::str_replace(stringr::fixed("*"), "") %>%
       stringr::str_replace(stringr::fixed(":"), "")]

       write_nuggets_nmers(dt, alleles)


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

    dtfn <-
      {
        data.table::rbindlist(
          list(
            dt[!netMHC %>% is.na &
            nmer_l < 15,
              .SD, .SDcols = c("netMHC", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHC, nmer_l) %>%
            write_netmhc_nmers("netMHC"),

            dt[!netMHCpan %>% is.na &
            nmer_l < 15,
              .SD, .SDcols = c("netMHCpan", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHCpan, nmer_l) %>%
            write_netmhc_nmers("netMHCpan"),

            dt[!netMHCII %>% is.na &
            nmer_l == 15,
              .SD, .SDcols = c("netMHCII", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHCII, nmer_l) %>%
            write_netmhc_nmers("netMHCII"),

            dt[!netMHCIIpan %>% is.na &
            nmer_l == 15,
              .SD, .SDcols = c("netMHCIIpan", "nmer", "nmer_l")] %>%
            data.table::copy %>%
            data.table::setkey(netMHCIIpan, nmer_l) %>%
            write_netmhc_nmers("netMHCIIpan")
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
              any){
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
            dt %>%
              data.table::setnames("allele", ptype)

            return(dt)
                     })
   return(dtl)
 }



## ---- write_nuggets_nmers
#' Internal function to output nmers for mhcnuggets prediction to disk
#'
#' @param dt Data table of nmers.
#' @param alleles Data table of 2 columns, 1. alleles properly formatted mhcnuggets.
#'
#' @export write_nuggets_nmers
#' @md

write_nuggets_nmers <- function(dt, alleles){

  if (dt %>% nrow == 0) return(NULL)

  if (!c("mhcnuggets", "nmer", "nmer_l") %chin% (dt %>% names) %>% any)
          stop("dt must contain mhcnuggets and nmer columns")


  # write tables for gru input

    if (dt[mhcnuggets %chin% alleles[type == "mhcnuggets_gru", allele]] %>% nrow > 0){

      mnug_dt <- dt[mhcnuggets %chin% alleles[type == "mhcnuggets_gru", allele] &
                    nmer_l <= 11,
                  .SD, .SDcols = c("mhcnuggets", "nmer")] %>%
               data.table::setnames(c("mhcnuggets", "nmer"),
                                    c("allele", "peptide")) %>%
               unique

      suppressWarnings(
           for (i in mnug_dt[, allele %>% unique]){
            breakpoints <- ((mnug_dt[allele == i] %>% nrow)/100) %>% ceiling

      parallel::mclapply(mnug_dt[allele == i] %>% split(1:breakpoints), function(dt){

            filename <- paste0("mhcnuggets_input_gru_",
                               i,
                               "_",
                               uuid::UUIDgenerate() %>%
                               substr(1, 18), ".csv")

            data.table::fwrite(dt[allele == i, peptide] %>%
                        data.table::as.data.table,
                        filename,
                        col.names = FALSE)
                        })
                      }
                  )
              }

  # write tables for lstm input

    if (dt[mhcnuggets %chin% alleles[type == "mhcnuggets_lstm", allele]] %>% nrow > 0){

      mnug_dt <- dt[mhcnuggets %chin% alleles[type == "mhcnuggets_lstm", allele] &
                    nmer_l <= 11,
                  .SD, .SDcols = c("mhcnuggets", "nmer")] %>%
        data.table::setnames(c("mhcnuggets", "nmer"),
                             c("allele", "peptide")) %>%
        unique

      suppressWarnings(
           for (i in mnug_dt[, allele %>% unique]){
            breakpoints <- ((mnug_dt[allele == i] %>% nrow)/100) %>% ceiling

       parallel::mclapply(mnug_dt[allele == i] %>% split(1:breakpoints), function(dt){

            filename <- paste0("mhcnuggets_input_lstm_",
                               i,
                               "_",
                               uuid::UUIDgenerate() %>%
                               substr(1, 18), ".csv")

            data.table::fwrite(dt[allele == i, peptide] %>%
                        data.table::as.data.table,
                        filename,
                        col.names = FALSE)
                        })
                      }
                  )
              }
  }



## ---- write_netmhc_nmers
#' Internal function to output nmers for netMHC prediction to disk
#'
#' @param dt Data table of nmers.
#' @param type Character vector. Name of program to format for.
#'
#' @export write_netmhc_nmers
#' @md

write_netmhc_nmers <- function(dt, type){
    if (dt %>% nrow == 0) return(NULL)
    if (!c("nmer", "nmer_l") %chin% (dt %>% names) %>% any)
          stop("dt must contain nmer and nmer_l columns")

    combs <- data.table::CJ(dt[, get(type)] %>% unique,
                            dt[, nmer_l] %>% unique)

  dto <- parallel::mclapply(1:nrow(combs), function(i){

      dts <- dt[get(type) == combs$V1[i] & nmer_l == combs$V2[i]] %>%
      unique

      # parallelize over 100 peptide chunks
      chunks <- ((dts %>% nrow)/100) %>% ceiling

  suppressWarnings(
  dto <- parallel::mclapply(dts %>% split(1:chunks), function(dtw){

        filename <- paste0(type, "_",
                    uuid::UUIDgenerate() %>% substr(1, 18), ".csv")


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
    )
      return(dto)

    }) %>% data.table::rbindlist

  return(dto)
  }



## ---- get_ss_str
#' Parallelized function to create a space-separated character string between two integers
#'
#' @param x Integer. Starting integer.
#' @param y Integer. Ending integer.
#'
#' @export get_ss_str
#' @md

get_ss_str <- function(x, y){

mcMap(function(x, y) (x %>% as.integer):(y %>% as.integer) %>%
          paste(collapse = " "), x, y) %>% unlist
    }



## ---- garnish_predictions
#' Perform neoepitope prediction.
#'
#' Perform ensemble neoepitope prediction on a data table of missense mutations, insertions, deletions or gene fusions using netMHC, mhcflurry, and mhcnuggets.
#' @param path Path to input table ([acceptable formats](https://cran.r-project.org/web/packages/rio/vignettes/rio.html#supported_file_formats)).
#' @param dt Data table. Input data table from `garnish_variants` or `garnish_jaffa`, or a data table in one of these forms:
#'
#'dt with transcript id:
#'
#'
#'     Column name                 Example input
#'
#'     sample_id                   sample_1
#'     ensembl_transcript_id       ENST00000311936
#'     cDNA_change                 c.718T>A
#'     MHC                         HLA-A02:01 HLA-A03:01
#'                                 H-2-Kb H-2-Kb
#'                                 HLA-DRB111:07 [second type]
#'
#'
#'dt with peptide:
#'
#'     Column name                 Example input
#'
#'     sample_id                   <same as above>
#'     pep_mut                     MTEYKLVVVDAGGVGKSALTIQLIQNHFV
#'     mutant_index                all
#'                                 7
#'                                 7 13 14
#'     MHC                         <same as above>
#'
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
#' * **pep_type**: type of peptide
#' * **pep_mut**: mutant peptide sequence
#' * **pep_wt**: wt peptide sequence
#' * **mismatch_s**: starting index of mutant peptide sequence
#' * **mismatch_l**: ending index of mutant peptide sequence
#' * **mutant_index**: index of mutant peptide
#' * **nmer**: nmer for prediction
#' * **nmer_i**: index of nmer in sliding window
#' * **\*_net**: netMHC prediction tool output
#' * **mhcflurry_\***: mhcflurry_ prediction tool output
#' * **mhcnuggets_\***: mhcnuggets_ prediction tool output
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
#'# run full pipeline
#'
#'    # download an example VCF
#'    dt <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'    # extract variants
#'    antigen.garnish::garnish_variants %>%
#'
#'    # add test MHC types
#'        .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
#'                    "H-2-Kb H-2-IAd",
#'                    "HLA-A*01:47 HLA-DRB1*03:08")] %>%
#'
#'    # predict neoepitopes
#'    antigen.garnish::garnish_predictions %>%
#'
#'    # summarize predictions
#'    antigen.garnish::garnish_summary %T>%
#'    print
#'
#'}
#'
#'\dontrun{
#'# input a data table of transcripts
#'
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
#'  antigen.garnish::garnish_predictions %T>%
#'  str
#' }
#'
#'\dontrun{
#'# input a data table of peptides for all MHC types
#'
#'library(data.table)
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  dt <- data.table::data.table(
#'           sample_id = "test",
#'           pep_mut = "MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP",
#'           mutant_index = "12",
#'           MHC = "all") %>%
#'  antigen.garnish::garnish_predictions %T>%
#'  str
#' }
#'
#'\dontrun{
#'# input from Microsoft Excel
#'
#'    # download an example Excel file
#'      path <- "antigen.garnish_test_input.xlsx" %T>%
#'        utils::download.file("http://get.rech.io/antigen.garnish_test_input.xlsx", .)
#'
#'    # predict neoepitopes
#'      dt <- antigen.garnish::garnish_predictions(path = path) %T>%
#'      str
#' }
#' @export garnish_predictions
#' @md

garnish_predictions <- function(dt = NULL,
                                path = NULL,
                               assemble = TRUE,
                               generate = TRUE,
                               predict = TRUE,
                               humandb = "GRCh38",
                               mousedb = "GRCm38"){

  # remove temporary files on exit
  on.exit({
    message("Removing temporary files")
    list.files(pattern = "(netMHC|mhcflurry|mhcnuggets).*_[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}\\.csv") %>% file.remove
  })

  if (missing(dt) & missing(path)) stop("dt and path are missing.")
  if (!missing(dt) & !missing(path)) stop("Choose dt or path input.")

  if (missing(dt) & !missing(path))
    dt <- rio::import(path) %>%
    data.table::as.data.table

  # specify transcript vs. direct cDNA / mutant index input
    if (c("sample_id", "ensembl_transcript_id", "cDNA_change", "MHC") %chin%
            (dt %>% names) %>% all) input_type <- "transcript"

    if (c("sample_id", "pep_mut", "MHC") %chin%
            (dt %>% names) %>% all) input_type <- "peptide"
            dt[, frameshift := FALSE]


  if (!exists("input_type"))
      stop("
Input data table must be from
garnish.variants or in
one of these forms:

dt with transcript id:


     Column name                 Example input

     sample_id                   sample_1
     ensembl_transcript_id       ENST00000311936
     cDNA_change                 c.718T>A
     MHC                         HLA-A02:01 HLA-A03:01
                                 H-2-Kb H-2-Kb
                                 HLA-DRB111:07 [second type]


dt with peptide:

     Column name                 Example input

     sample_id                   <same as above>
     pep_mut                     MTEYKLVVVDAGGVGKSALTIQLIQNHFV
     mutant_index                all
                                 7
                                 7 13 14
     MHC                         <same as above>
      ")

if (assemble & input_type == "transcript"){

    dt %<>% get_metadata(humandb = humandb, mousedb = mousedb)

    # extract cDNA changes and make cDNA
    dt %<>% extract_cDNA
    dt %<>% make_cDNA

    # translate protein sequences
    dt[, pep_wt := coding %>% translate_cDNA]
    dt[, pep_mut := coding_mut %>% translate_cDNA]
    dt[, frameshift := FALSE]

    dt[, cDNA_delta := ((coding_mut %>% nchar) - (coding %>% nchar)) / 3 ]

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

    # create a space-separated vector of mutant indices
        dt[, mutant_index := mismatch_s %>% as.character]
        dt[mismatch_l > mismatch_s, mutant_index :=
            get_ss_str(mismatch_s, mismatch_l)]

    # warning if mutant_index == NA
      if (dt[is.na(mutant_index) | (nchar(pep_mut) - as.numeric(mutant_index) < 1)] %>%
            nrow > 1){
        failn <- dt[is.na(mutant_index) | (nchar(pep_mut) - as.numeric(mutant_index) < 1)] %>% nrow
        warning(paste0("Could not determine mutant index for ", failn, " records."))
      }
  }

if (assemble & input_type == "peptide"){
    dt[mutant_index == "all", mutant_index :=
      get_ss_str(1, pep_mut %>% nchar)]
}

if (generate){

  message("Generating variants")

  # generation a uuid for each unique variant

   suppressWarnings(dt[, var_uuid :=
      parallel::mclapply(1:nrow(dt),
      uuid::UUIDgenerate) %>% unlist])

  # separate over mutant indices

    if (dt[, mutant_index %>% unique] %>%
        stringr::str_detect(" ") %>% any){
      dts <- dt %>% tidyr::separate_rows("mutant_index", sep = " ")
    } else {
      dts <- dt
    }

  # convert back to numeric

    dts[, mutant_index := mutant_index %>% as.numeric]

  # generate a data table of unique variants for peptide generation

    dtnfs <- dts[frameshift == FALSE]
    dtfs <- dts[frameshift == TRUE]

   if (input_type == "transcript"){
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
            .[, pep_type := "mut_other"])) %>%

         # take unique peptides
            .[, .SD, .SDcols = c("var_uuid",
                                 "pep_type",
                                 "pep_base",
                                 "mutant_index")] %>%
            unique
    }

   if (input_type == "peptide"){
      basepep_dt <- dtnfs %>%
                      data.table::copy %>%
                      .[, pep_base := pep_mut] %>%
                      .[, pep_type := "mut_other"]
     }

  # filter mutant_index
    basepep_dt <- basepep_dt[(nchar(pep_base) - mutant_index) >= 1]

  if (basepep_dt %>% nrow == 0)
    return("no variants for peptide generation")

    sink(file = "/dev/null")
    nmer_dt <- make_nmers(basepep_dt) %>% .[, nmer_l := nmer %>% nchar]
    sink()
    dt <- merge(dt, nmer_dt,
            by = intersect(names(dt), names(nmer_dt)),
            all.x = TRUE)

    ## drop out single wt nmer from rolling window over fusion peptides from JAFFA input
    if("fus_tx" %chin% names(dt)) dt <- dt %>%
      .[, drop := stringi::stri_detect_fixed(pattern = nmer, str = pep_gene_1)] %>%
                  .[drop == FALSE] %>%
                    .[, drop := NULL]


    # generation a uuid for each unique nmer
    suppressWarnings(dt[, nmer_uuid :=
                    parallel::mclapply(1:nrow(dt),
                    uuid::UUIDgenerate) %>% unlist])
    if (input_type == "transcript"){

         dt %<>% make_DAI_uuid

        # remove mut == wt by sample_id
        # remove wt missing mut counterpart

        for (id in dt[, sample_id %>% unique]){

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

    # DAI cannot be calculated with peptide input
    if (input_type == "peptide")
      dt[, dai_uuid := NA %>% as.character]

}

if (predict){

  dtl <- dt %>% get_pred_commands

  # run commands

   if (check_pred_tools() %>% unlist %>% all){

      dto <- run_netMHC(dtl[[2]])

      run_mhcflurry()
      run_mhcnuggets()

      dt <- merge_predictions(dto, dtl[[1]])

      cols <- dt %>% names %include% "(best_netMHC)|(mhcflurry_prediction$)|(mhcnuggets_pred_gru)|(mhcnuggets_pred_lstm)"

      confi <- function(dt){

        dtl <- parallel::mclapply(1:nrow(dt), function(i){

          if (dt[i ,] %>% unlist %>% na.omit %>% unique %>% length <= 1){
            return(list(NA, NA))}
            t.test(dt[i ,])$conf.int[1:2] %>% as.list
                          })

        lo <- lapply(dtl, function(x){x[[1]]}) %>% unlist

        up <- lapply(dtl, function(x){x[[2]]}) %>% unlist
        return(list(lo, up))
                            }

      dt[, c("Lower.CI", "Upper.CI") := confi(.SD), .SDcols = cols]
      } else {
        warning("Missing prediction tools in PATH, returning without predictions.")
      }

   }
   return(dt)
}



## ---- make_nmers
#' Internal function for parallelized `nmer` creation.
#'
#' @param dt Data table. Input data table from `garnish_predictions`.
#'
#' @export make_nmers
#' @md

make_nmers <- function(dt){


  if (!c("var_uuid",
         "pep_type",
         "pep_base",
         "mutant_index") %chin% (dt %>% names) %>% all)
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

      nmer_dt <- lapply((15:8), function(pl){

        mut_frag_t <- dt$pep_base[n] %>% strsplit("",
                            fixed = TRUE) %>% unlist
        mut_frag_loc <- dt$mutant_index[n]

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
