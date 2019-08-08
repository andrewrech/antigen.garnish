


#' Internal function to Smith-Waterman align two vectors of peptides.
#'
#' @param col1 Character. Vector of sequences to align to col2.
#' @param col2 Character. Vector of sequences to align to col1.
#' @param gap_open Numeric. Cost to open gapped alignment. Default is -11.
#' @param gap_extend Numeric. Cost to extend a gap. Default is -1.
#'
#' @export SW_align
#' @md
SW_align <- function(col1,
                      col2,
                      gap_open = -11,
                      gap_extend = -1){

        al <- Biostrings::pairwiseAlignment(col1, col2,
                                  substitutionMatrix = "BLOSUM62",
                                  gapOpening = gap_open,
                                  gapExtension = gap_extend,
                                  type = "local",
                                  scoreOnly = TRUE)

        if (length(al) == 0) al <- as.numeric(NA)

        return(al)

}

#' Internal function to run run_mhcflurry commands.
#'
#' @export run_mhcflurry
#' @md

run_mhcflurry <- function(){

  message("Running mhcflurry in parallel")

    list.files(pattern = "mhcflurry_input.*csv") %>%

mclapply(., function(x){
      paste0("mhcflurry-predict ", x, " --out ", x %>%
            stringr::str_replace("input", "output")) %>%
      system
          })
}




#' Internal function to run run_mhcnuggets commands.
#'
#' @export run_mhcnuggets
#' @md

run_mhcnuggets <- function(){

  gruf <- list.files(pattern = "mhcnuggets_input_gru.*csv")

  if (gruf %>% length > 1)
    message("Running mhcnuggets with -m gru")

    pypath  <- system("which predict.py", intern = TRUE)

parallel::mclapply(gruf, function(x){
        paste0("python ", pypath, " -m gru -w ", dirname(pypath), "/../saves/production/mhcnuggets/",
               stringr::str_extract(string = x, pattern = "(?<=_)H.*(?=_)"), ".h5 -p ", x,
               " > ",
               x %>% stringr::str_replace("input", "output")) %>%
          system

      return(NULL)

      })



  lf <- list.files(pattern = "mhcnuggets_input_lstm.*csv")

  if (lf %>% length > 1)
  message("Running mhcnuggets with -m lstm")

  pypath  <- system("which predict.py", intern = TRUE)

parallel::mclapply(lf, function(x){
        paste0("python ", pypath, " -m lstm -w ",  dirname(pypath), "/../saves/production/mhcnuggets_beta/",
               stringr::str_extract(string = x, pattern = "(?<=_)H.*(?=_)"), ".h5 -p ", x,
               " > ",
               x %>% stringr::str_replace("input", "output")) %>%
          system

  return(NULL)

      })
}




#' Internal function to run netMHC commands.
#'
#' @param dt Data table of prediction commands to run.
#'
#' @export run_netMHC
#' @md

run_netMHC <- function(dt){

  if (!"command" %chin% (dt %>% names))
    stop("dt must contain command column")

  message("Running netMHC in parallel")

  fn <- paste(Sys.time() %>% stringr::str_replace_all("\\ |-|:", ""),
  "_netMHC_commands.txt", sep = "")

  # set syntax to avoid invalid internal selfref warning
  set(dt,
    j = "outname",
    value = dt[, filename] %>%
    stringr::str_replace("\\.csv$", "_o.csv"))


  dt[, command2 := paste(command,
                              " > ",
                              outname,
                              sep = "")]

  dt[, .SD, .SDcols = "command2"] %>%
    data.table::fwrite(fn, col.names = FALSE)

  system(paste("parallel --bar --no-run-if-empty --joblog parallel_log.txt <",
    fn, sep = " "))

  outfiles <- dt[, outname]

  # check for expected output from parallel run
  if (!all(file.exists(outfiles)))
    stop("netMHC did not produce output for all commands.")

  logf <- data.table::fread("parallel_log.txt")

  if (any(logf[, Exitval != 0])){

    cmd <- logf[Exitval != 0, Command] %>% paste(collapse  =  "\n")

    stop(paste("netMHC errored on a command, check MHC syntax with list_mhc\\(\\)",
    "and make sure netMHC is in path. Errored command was:",
    cmd, sep = " "))

  }

  file.remove(fn)

  file.remove("parallel_log.txt")

  # generate command/output file list
  esl <- lapply(1:nrow(dt), function(n){

    return(list(dt[n][, command], dt[n][, outname]))

  })

  dtl <- esl %>% collate_netMHC

      return(dtl)

}


#' Print data directory error
#'
#' @export ag_data_err
#' @md

ag_data_err <- function(){
    err <- paste(
    "",
    "Unable to locate antigen.garnish data directory,",
    "or directory contents are missing.",
    "",
    "Paths searched are $AG_DATA_DIR, $HOME",
    "and the current working directory.",
    "",
    "To set a custom path to the antigen.garnish data folder",
    "set environomental variable AG_DATA_DIR from the shell",
    "or from R using Sys.setenv",
    "",
    "Re-download installation data:",
    '$ curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz',
    "",
    "Documentation:",
    "https://neoantigens.io",
    sep = "\n"
    )

  stop(err)

}

#' Internal function to check for netMHC tools, mhcflurry and mhcnuggets
#'
#' @param ag_dirs Character vector. Directories to check for `antigen.garnish` data directory.
#'
#' @export check_pred_tools
#' @md

check_pred_tools <- function(ag_dirs = c(
                                          paste0(
                                            getwd(),
                                            "/",
                                            "antigen.garnish"),
                                          paste0(
                                            Sys.getenv("HOME"),
                                            "/",
                                            "antigen.garnish")
                                            )){


  # vector of directories to look in for data

 if (!Sys.getenv("AG_DATA_DIR") == "")
   if (!Sys.getenv("AG_DATA_DIR") %>% dir.exists){

     message("$AG_DATA_DIR does not exist")
     ag_data_err()

   }

   if (Sys.getenv("AG_DATA_DIR") %>% dir.exists){

     ag_dir <- Sys.getenv("AG_DATA_DIR")

   }


 if (Sys.getenv("AG_DATA_DIR") == ""){

    for(i in ag_dirs)
      if (i %>% dir.exists){
        ag_dir <- i
        break
      }

    Sys.setenv(AG_DATA_DIR = ag_dir)

   }

  default_path <- paste0(ag_dir,
                  c(
                    "/netMHC/netMHC-4.0",
                    "/netMHC/netMHCII-2.2",
                    "/netMHC/netMHCIIpan-3.1",
                    "/netMHC/netMHCpan-3.0",
                    "/mhcnuggets/scripts"
                    )) %>%
                      normalizePath %>%
                      paste(collapse = ":")

  tool_status <- list(
                mhcflurry   = TRUE,
                netMHC      = TRUE,
                netMHCpan   = TRUE,
                netMHCII    = TRUE,
                netMHCIIpan = TRUE,
                mhcnuggets  = TRUE
                )

  # conditional to prevent infinitely growing path when running in long loops
  if(!grepl(pattern = default_path, x = Sys.getenv("PATH"), fixed = TRUE))
    Sys.setenv(PATH = paste0(default_path, ":", Sys.getenv("PATH")))

    if (suppressWarnings(system('which mhcflurry-predict 2> /dev/null', intern = TRUE)) %>%
          length == 0){
          message("mhcflurry-predict is not in PATH\n       Download: https://github.com/hammerlab/mhcflurry")
        tool_status$mhcflurry <- FALSE
        }
    if (suppressWarnings(system('which netMHC 2> /dev/null', intern = TRUE)) %>%
          length == 0){
            message("netMHC is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHC/")
          tool_status$netMHC <- FALSE
          }
    if (suppressWarnings(system('which netMHCpan 2> /dev/null', intern = TRUE)) %>%
          length == 0){
            message("netMHCpan is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCpan/")
          tool_status$netMHCpan <- FALSE
          }
    if (suppressWarnings(system('which netMHCII 2> /dev/null', intern = TRUE)) %>%
          length == 0){
            message("netMHCII is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCII/")
          tool_status$netMHCII <- FALSE
          }
    if (suppressWarnings(system('which netMHCIIpan 2> /dev/null', intern = TRUE)) %>%
          length == 0){
            message("netMHCIIpan is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCIIpan/")
          tool_status$netMHCIIpan <- FALSE
         }
    if (suppressWarnings(system('which predict.py 2> /dev/null', intern = TRUE)) %>%
          length == 0){
            message("mhcnuggets predict.py is not in PATH\n       Download: https://github.com/KarchinLab/mhcnuggets")
          tool_status$mhcnuggets <- FALSE
         }

          return(tool_status)
}
