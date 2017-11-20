


## ---- run_mhcflurry
#' Internal function to run run_mhcflurry commands.
#'
#' @export run_mhcflurry
#' @md

run_mhcflurry <- function(){

  message("Running mhcflurry in parallel")

    list.files(pattern = "mhcflurry_input.*csv") %>%

  mclapply(., function(x){
      paste0("mhcflurry-predict ", x, " > ", x %>%
            stringr::str_replace("input", "output")) %>%
      system
          })
}



## ---- run_mhcnuggets
#' Internal function to run run_mhcnuggets commands.
#'
#' @export run_mhcnuggets
#' @md

run_mhcnuggets <- function(){

  gruf <- list.files(pattern = "mhcnuggets_input_gru.*csv")

  if (gruf %>% length > 1)
    message("Running mhcnuggets with -m gru")

  parallel::mclapply(gruf, function(x){
        paste0("python $(which predict.py) -m gru -w $(dirname $(which predict.py))/../saves/production/mhcnuggets/",
               stringr::str_extract(string = x, pattern = "(?<=_)H.*(?=_)"), ".h5 -p ", x,
               " > ",
               x %>% stringr::str_replace("input", "output")) %>%
          system
      })



  lf <- list.files(pattern = "mhcnuggets_input_lstm.*csv")

  if (lf %>% length > 1)
  message("Running mhcnuggets with -m lstm")

  parallel::mclapply(lf, function(x){
        paste0("python $(which predict.py) -m lstm -w $(dirname $(which predict.py))/../saves/production/mhcnuggets_beta/",
               stringr::str_extract(string = x, pattern = "(?<=_)H.*(?=_)"), ".h5 -p ", x,
               " > ",
               x %>% stringr::str_replace("input", "output")) %>%
          system
      })
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

  message("Running netMHC in parallel")

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



## ---- check_pred_tools
#' Internal function to check for netMHC tools, mhcflurry and mhcnuggets
#'
#' @export check_pred_tools
#' @md

check_pred_tools <- function(){

  default_path <- paste0(system('echo /usr/bin/local/', intern = TRUE),
                  c(
                  "/antigen.garnish/netMHC/netMHC-4.0",
                  "/antigen.garnish/netMHC/netMHCII-2.2",
                  "/antigen.garnish/netMHC/netMHCIIpan-3.1",
                  "/antigen.garnish/netMHC/netMHCpan-3.0",
                  "/antigen.garnish/mhcnuggets/scripts"
                  )) %>%
                  paste(collapse = ":")

  tool_status <- list(
                mhcflurry = TRUE,
                netMHC = TRUE,
                netMHCpan = TRUE,
                netMHCII = TRUE,
                netMHCIIpan = TRUE,
                mhcnuggets = TRUE
                )

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
