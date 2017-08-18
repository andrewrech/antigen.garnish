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
        length == 0){
        message("mhcflurry-predict is not in PATH\n       Download: https://github.com/hammerlab/mhcflurry")
      PATH_status$mhcflurry <- FALSE
      }
  if (suppressWarnings(system('which netMHC 2> /dev/null', intern = TRUE)) %>%
        length == 0){
          message("netMHC is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHC/")
        PATH_status$netMHC <- FALSE
        }
  if (suppressWarnings(system('which netMHCpan 2> /dev/null', intern = TRUE)) %>%
        length == 0){
          message("netMHCpan is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCpan/")
        PATH_status$netMHCpan <- FALSE
        }
  if (suppressWarnings(system('which netMHCII 2> /dev/null', intern = TRUE)) %>%
        length == 0){
          message("netMHCII is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCII/")
        PATH_status$netMHCII <- FALSE
        }
  if (suppressWarnings(system('which netMHCIIpan 2> /dev/null', intern = TRUE)) %>%
        length == 0){
          message("netMHCIIpan is not in PATH\n       Download: http://www.cbs.dtu.dk/services/NetMHCIIpan/")
        PATH_status$netMHCIIpan <- FALSE
        }
        return(PATH_status)
}
