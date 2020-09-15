


#' Internal function to Smith-Waterman align two vectors of peptides.
#'
#' @param col1 Character. Vector of sequences to align to col2.
#' @param col2 Character. Vector of sequences to align to col1.
#' @param gap_open Numeric. Cost to open gapped alignment. Default is -11.
#' @param gap_extend Numeric. Cost to extend a gap. Default is -1.
#'
#' @md
make_sw_alignment <- function(col1,
                              col2,
                              gap_open = -11,
                              gap_extend = -1) {
  al <- Biostrings::pairwiseAlignment(col1, col2,
    substitutionMatrix = "BLOSUM62",
    gapOpening = gap_open,
    gapExtension = gap_extend,
    type = "local",
    scoreOnly = TRUE
  )

  if (length(al) == 0) al <- as.numeric(NA)

  return(al)
}

#' Internal function to run run_mhcflurry commands.
#'
#' @md

run_mhcflurry <- function() {
  message("Running mhcflurry in parallel")

  fn <- list.files(pattern = "mhcflurry_input.*csv")

  fn %>%
    seq_along() %>%
    mclapply(function(x) {
      message(paste0("Input file ", x, " of ", length(fn)))

      paste0("mhcflurry-predict ", fn[x], " --out ", fn[x] %>%
        stringr::str_replace("input", "output"), " 2>>mhcflurry-predict.log") %>%
        system()
    })
}

#' Internal function to run netMHC commands.
#'
#' @param dt Data table of prediction commands to run.
#'
#' @md

run_netMHC <- function(dt) {
  if (!"command" %chin% (dt %>% names())) {
    stop("dt must contain command column")
  }

  fn <- paste(Sys.time() %>% stringr::str_replace_all("\\ |-|:", ""),
    "_netMHC_commands.txt",
    sep = ""
  )

  scriptn <- fn %>% stringr::str_replace("\\.txt$", ".sh")

  plog <- fn %>% stringr::str_replace("\\.txt$", "_parallel.log")

  # set syntax to avoid invalid internal selfref warning
  set(dt,
    j = "outname",
    value = dt[, filename] %>%
      stringr::str_replace("\\.csv$", "_o.csv")
  )


  dt[, command2 := paste(command,
    " > ",
    outname,
    sep = ""
  )]

  # write 10 commands per line sep by ;
  cmdlist <- dt[, .SD %>% unique(), .SDcols = "command2"] %>% split(1:(nrow(.) / 10))

  cmdlist <- lapply(cmdlist, function(i) {
    i %>%
      unlist() %>%
      paste(collapse = ";")
  }) %>%
    unlist() %>%
    unique() %>%
    as.data.table()

  cmdlist %>% data.table::fwrite(fn, col.names = FALSE)

  message("Running netMHC in parallel")

  # we write to file via R base cat fxn (echo is bash only)
  # set jobs max value 90% of available CPUs, delay 0.2 on spawning jobs for safety
  # 5000 is a conservative line length to be safe, could go much higher  probably
  # memory is not protected here, user will need to chunk if RAM limited
  cat(
    paste("parallel --eta -s 5000 --delay 0.1 --no-run-if-empty --joblog ",
      plog, " < ", fn,
      sep = ""
    ),
    file = scriptn,
    fill = TRUE
  )

  Sys.chmod(scriptn, mode = "0777")

  # try to run with sh, if sh not present, will let system run default shell
  shpath <- system("which sh", intern = TRUE)

  # replace whitespace if shpath is empty to let default shell run it
  cmd <- paste(shpath, scriptn, sep = " ") %>%
    stringr::str_replace("^\\ ", "./")

  # parallel will run commands with whatever shell script is invoked with
  exstat <- system(cmd, wait = TRUE) %>% unlist()

  outfiles <- dt[, outname]

  logf <- data.table::fread(plog)

  # check all commands exited with 0 from parallel logfile
  if (any(logf[, Exitval != 0])) {
    cmd <- logf[Exitval != 0, Command] %>% paste(collapse = "\n")

    stop(paste("netMHC errored on a command, check MHC syntax with `list_mhc()`",
      "and make sure netMHC is in path. Errored command was:",
      cmd,
      sep = " "
    ))
  }

  # check for expected output files from parallel run
  if (!all(file.exists(outfiles))) {
    stop("netMHC did not produce output for all commands.")
  }

  # just in case files get created and parallel logfile doesnt catch error, check here.
  if (exstat != 0) stop("Exit status from running netMHC parallel script is not zero.")

  file.remove(fn)

  file.remove(scriptn)

  file.remove(plog)

  # generate command/output file list
  esl <- lapply(1:nrow(dt), function(n) {
    return(list(dt[n][, command], dt[n][, outname]))
  })

  dtl <- esl %>% collate_netMHC()

  return(dtl)
}

#' Print data directory error
#'
#' @md


.ag_data_err <- function() {
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
#' @md

configure_netMHC_tools <- function(dir) {
  owd <- getwd()

  on.exit({
    setwd(owd)
  })

  # netMHC parent dir in ag data dir
  npd <- file.path(dir, "netMHC")

  if (!dir.exists(npd)) {
    stop("netMHC parent directory in antigen.garnish data directory not found.")
  }

  # get path to scripts
  setwd(npd)

  tool_paths <- list.files()

  f <- lapply(tool_paths, function(i) {

    # pattern to find netMHC scripts
    list.files(i,
      pattern = "netMHC(II)?(pan)?-?([0-9]\\.[0-9])?$",
      full.names = TRUE, ignore.case = TRUE
    )
  }) %>% unlist()

  message("Checking netMHC scripts in antigen.garnish data directory.")
  # sed scripts to correct paths
  # install data from DTU
  io <- lapply(f, function(i) {

    # rename to take off version, necessary because commands are built into table
    # in get_pred_commands and check_pred_tools hasn't run at that point
    io <- file.path(
      dirname(i),
      basename(i) %>% stringr::str_replace("-.*$", "")
    )

    line <- file.path(dir, "netMHC", dirname(i))

    # downloading the data takes a long time, must check if loop already ran
    catfile <- system(paste("cat", i), intern = TRUE) %>%
      stringr::str_replace_all("/", "")

    # adjust patterns to have no slashes
    p1 <- line %>% stringr::str_replace_all("/", "")

    # if correct directory is already in the script, skip
    if (any(grepl(catfile, pattern = p1))) {
      return(io)
    }

    # properly escape for sed call
    line <- paste("setenv NMHOME ", line %>% stringr::str_replace_all("/", "\\\\/"), sep = "")
    line2 <- paste("setenv TMPDIR ", "$HOME", sep = "")

    # formatting for these links is  not consistent,  try no uname first
    # netMHC is always capitalized in the links it seems
    # dirname(i) will always have version number
    link_to_data <- paste0(
      "http://www.cbs.dtu.dk/services/",
      dirname(i) %>% stringr::str_replace("net", "Net"),
      "/data.tar.gz"
    )

    cmd <- paste("curl -fsSL", link_to_data, ">", "dtu.tar.gz")

    suppressWarnings({
      curl_status <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    })

    if (!is.null(attr(curl_status, "status")) &&
      attr(curl_status, "status") == 22) {
      un <- system("uname -s", intern = TRUE)

      if (!un %chin% c("Linux")) stop("Linux only.")

      cmd <- cmd %>%
        stringr::str_replace("data\\.tar\\.gz", paste0("data.", un, ".tar.gz"))

      message(paste("Downloading data for", dirname(i), "on", un))

      suppressWarnings({
        curl_status <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
      })
    }

    if (!file.exists("dtu.tar.gz") || file.info("dtu.tar.gz")$size == 0) {
      stop(paste("Unable to download data tar from DTU. See ReadMe in", file.path(dir, "netMHC", dirname(i))))
    }

    # move it into the right folder and untar
    dtar <- file.path(dirname(i), "dtu.tar.gz")

    file.rename("dtu.tar.gz", dtar)

    setwd(dirname(i))

    system(paste("tar -xzvf", basename(dtar)))

    setwd(npd)

    # versionless backup file
    bk <- paste0(io, ".bk")

    if (!file.exists(bk)) file.rename(i, bk)

    cmd <- paste("cat ", bk, " | ",
      # place holders for variable references
      "sed 's/\\$NMHOME/placeholderxxx/' | ",
      "sed 's/\\$TMPDIR/placeholderyyy/' | ",
      "sed 's/setenv.*NMHOME.*$/", line, "/' | ",
      # reset placeholder
      "sed 's/placeholderxxx/$NMHOME/' | ",
      "sed 's/setenv.*TMPDIR.*$/", line2, "/' | ",
      # reset placeholder
      "sed 's/placeholderyyy/$TMPDIR/' > ", io,
      sep = ""
    )

    system(cmd)

    system(paste("chmod u+x", io))

    return(io)
  }) %>% unlist()

  setwd(owd)

  return(io)
}

#' Internal function to check for netMHC tools and mhcflurry
#'
#' @param ag_dirs Character vector. Directories to check for `antigen.garnish` data directory.
#'
#' @md

check_pred_tools <- function(ag_dirs = c(
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
                             )) {


  # vector of directories to look in for data

  if (!Sys.getenv("AG_DATA_DIR") == "") {
    message("Environmental variable AG_DATA_DIR for the antigen.garnish data directory is unset. Checking standard directories.")
  }

  if (!Sys.getenv("AG_DATA_DIR") == "") {
    if (!Sys.getenv("AG_DATA_DIR") %>% dir.exists()) {
      message(Sys.getenv("AG_DATA_DIR"), " does not exist.")
      .ag_data_err()
    }
  }

  if (Sys.getenv("AG_DATA_DIR") %>% dir.exists()) {
    ag_dir <- Sys.getenv("AG_DATA_DIR")
  }

  if (Sys.getenv("AG_DATA_DIR") == "") {
    for (i in ag_dirs) {
      if (i %>% dir.exists()) {
        ag_dir <- i
        break
      }
    }

    Sys.setenv(AG_DATA_DIR = ag_dir)
  }

  scripts <- configure_netMHC_tools(ag_dir)

  default_path <- paste0(
    ag_dir,
    c(
      file.path(
        "/netMHC",
        list.files(file.path(ag_dir, "netMHC"))
      ),
      file.path("/ncbi-blast-2.10.1+-src/c++/ReleaseMT/bin")
    )
  ) %>%
    normalizePath() %>%
    paste(collapse = ":")

  tool_status <- list(
    mhcflurry   = TRUE,
    netMHC      = TRUE,
    netMHCpan   = TRUE,
    netMHCII    = TRUE,
    netMHCIIpan = TRUE,
    blastp      = TRUE
  )

  # conditional to prevent infinitely growing path when running in long loops
  if (!grepl(pattern = default_path, x = Sys.getenv("PATH"), fixed = TRUE)) {
    Sys.setenv(PATH = paste0(default_path, ":", Sys.getenv("PATH")))
  }

  if (suppressWarnings(system("which mhcflurry-predict 2> /dev/null", intern = TRUE)) %>%
    length() == 0) {
    message("mhcflurry-predict is not in PATH\n       Download: https://github.com/hammerlab/mhcflurry")
    tool_status$mhcflurry <- FALSE
  }

  if (suppressWarnings(system("which blastp 2> /dev/null", intern = TRUE)) %>%
    length() == 0) {
    message("blastp is not in PATH\n       Download: https://blast.ncbi.nlm.nih.gov/Blast.cgi")
    tool_status$blastp <- FALSE
  }

  lapply(scripts, function(i) {
    f <- basename(i)

    if (suppressWarnings(system(paste("which", f, "2> /dev/null"), intern = TRUE)) %>%
      length() == 0) {
      message(paste(f, " is not in PATH\n       Download: http://www.cbs.dtu.dk/services/"), sep = "")
      tool_status[[f]] <- FALSE
    }
  })

  return(tool_status)
}

#' Convenience inflix operator to return vector elements matching a regular expression.
#'
#' @param vector Vector.
#' @param pattern Pattern.
#'
#' @return A vector of elements in \code{vector} matching \code{pattern}.
#'
#' @export %include%
#' @md

`%include%` <- function(vector, pattern) {
  lv <- stringr::str_detect(vector, pattern)
  return(vector[lv])
}


#' Convenience inflix operator to return vector elements excluding those matching a regular expression.
#'
#' @param vector Vector.
#' @param pattern Pattern.
#'
#' @return A vector of elements in \code{vector} not matching \code{pattern}.
#'
#' @export %exclude%

`%exclude%` <- function(vector, pattern) {
  lv <- stringr::str_detect(vector, pattern)
  return(vector[!lv])
}
