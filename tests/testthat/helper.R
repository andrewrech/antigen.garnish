
# load antigen.garnish test helper functions

skip_pred_tools <- function() {

         if (!check_pred_tools() %>% unlist %>% all){
          testthat::skip("Skipping run_netMHC because prediction tools are not in PATH")
          }

         if (suppressWarnings(system('which blastp 2> /dev/null', intern = TRUE)) %>%
          length == 0)
          testthat::skip("Skipping garnish_affinity README example from VCF without fitness because ncbiblast+ is not in PATH")
        }

test_runner <- function(fn){
    cdir <- getwd()
    temp_dir <- uuid::UUIDgenerate()
    dir.create(temp_dir)
    setwd(temp_dir)
    fn()
    on.exit({
              setwd(cdir)
              unlink(paste0(cdir, "/", temp_dir),
                     force = TRUE,
                     recursive = TRUE)
              })
  }