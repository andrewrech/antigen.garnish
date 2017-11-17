#!/usr/bin/env sh

# install antigen.garnish R package and dependencies
# https://github.com/andrewrech/antigen.garnish

# install dependencies

  Rscript -e \
  'install.packages("devtools"); devtools::install_github("hadley/devtools"); devtools::install_github("r-lib/usethis"); install.packages("testthat")'

  Rscript -e \
  'source("https://bioconductor.org/biocLite.R"); biocLite(c("BiocInstaller", "ShortRead", "biomaRt", "ShortRead"), ask = FALSE)'

  Rscript -e \
  'remove.packages("data.table"); install.packages("data.table", type = "source", repos = "http://Rdatatable.github.io/data.table")'

  Rscript -e \
  'devtools::install_github(c("tidyverse/magrittr", "andrewrech/dt.inflix"))'