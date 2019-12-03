


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

  # run commands
  esl <- lapply(
         dt[, command],

function(command){
          # gen temp file and write to disk to save memory usage
          es <- uuid::UUIDgenerate() %>% stringr::str_replace_all("-", "")

          es <- paste(es, ".txt", sep = "")

          command2 <- paste(command, " > ", es)

          # run command
          system(command2)

          # if error, return empty dt
          if (!file.exists(es) || file.info(es)$size == 0)
              return(data.table::data.table(status = 1, command = command))

          return(list(command, es))

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

#' Internal function to configure netMHC suite of tools
#'
#' @param dir Character vector. Path to `antigen.garnish` data directory.
#'
#' @export configure_netMHC_tools
#' @md

configure_netMHC_tools <- function(dir){

  owd <- getwd()

  on.exit({setwd(owd)})

  # netMHC parent dir in ag data dir
  npd <- file.path(dir, "netMHC")

  if (!dir.exists(npd)) dir.create(npd)

  l <- list.files(npd,
  pattern = "netMHC.*tar")

  if (length(l) != 0){

    message("Untarring netMHC tools")

    owd <- getwd()

    setwd(npd)

    lapply(l, function(i){

      system(paste("tar -xzvf", i))

  })

  # remove so we don't do this with multiple runs
  # user was instructed to copy from home dir so originals are intact
  file.remove(l)

  setwd(owd)

  }

  # get path to scripts
  owd <- getwd()

  setwd(npd)

  tool_paths <- list.files()

  f <- lapply(tool_paths, function(i){

    # this pattern could one day break.
    list.files(i, pattern = "netMHC(II)?(pan)?-?([0-9]\\.[0-9])?$",
        full.names = TRUE)

  }) %>% unlist

  message("Changing permissions on all tools.")

  Sys.chmod(f, mode = "0777")

  message("Checking netMHC scripts in antigen.garnish data directory.")
  # sed scripts to correct paths
  io <- lapply(f, function(i){

    # check if scripts were already edited
    if (file.exists(file.path(dirname(i), "itwasedited.txt")))
      return(NULL)

    line <- file.path(dir, i)

    line <- paste("setenv NMHOME ", line %>% stringr::str_replace_all("/", "\\\\/"), sep = "")
    line2 <- paste("setenv TMPDIR ", owd %>% stringr::str_replace_all("/", "\\\\/"))

    cmd <- paste("cat ", i, " | ",
    # place holders for variable references
    "sed 's/\\$NMHOME/placeholderxxx/' | ",
    "sed 's/\\$TMPDIR/placeholderyyy/' | ",
    "sed 's/setenv.*NMHOME.*$/", line, "/' | ",
    # reset placeholder
    "sed 's/placeholderxxx/$NMHOME/' | ",
    "sed 's/setenv.*TMPDIR.*$/", line2, "/' | ",
    # reset placeholder
    "sed 's/placeholderyyy/$TMPDIR/' > ", "placeholder_script.txt", sep = "")

    system(cmd)

    # add marker that it was edited
    system(paste("touch", file.path(dirname(i), "itwasedited.txt")))

    file.remove(i)

    # rename to take off version, necessary because commands are built into table
    # in get_pred_commands and check_pred_tools hasn't run at that point
    io <- file.path(dirname(i),
    basename(i) %>% stringr::str_replace("-.*$", ""))

    file.rename("placeholder_script.txt", io)

    return(io)

  }) %>% unlist

  message("Done.")

  setwd(owd)

  return(io)

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

  scripts <- configure_netMHC_tools(ag_dir)

  default_path <- paste0(ag_dir,
                  c(file.path("/netMHC",
                  list.files(file.path(ag_dir, "netMHC"))),
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

    lapply(scripts, function(i){

      f <- basename(i)

      if (suppressWarnings(system(paste("which", f, "2> /dev/null"), intern = TRUE)) %>%
            length == 0){
              message(paste(f, " is not in PATH\n       Download: http://www.cbs.dtu.dk/services/"), sep = "")
            tool_status[[f]] <- FALSE
      }

    })

    if (suppressWarnings(system('which predict.py 2> /dev/null', intern = TRUE)) %>%
          length == 0){
            message("mhcnuggets predict.py is not in PATH\n       Download: https://github.com/KarchinLab/mhcnuggets")
          tool_status$mhcnuggets <- FALSE
         }

          return(tool_status)
}
