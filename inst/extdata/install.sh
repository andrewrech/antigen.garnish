#!/usr/bin/env sh

# install antigen.garnish R package and dependencies
# https://github.com/andrewrech/antigen.garnish

# check script tools

  OS=`uname`
  if [ ! "$OS" = "Linux" ]
  then
      echo "Linux only."
      exit 1
  fi

  if [ ! -x `which tar` ]; then
  echo "No suitable tar program found."
  exit 1
  fi

  if [ ! -x `which Rscript` ]; then
  echo "No suitable Rscript program found."
  exit 1
  fi

  if [ ! -x `which pip` ]; then
  echo "No suitable pip program found."
  exit 1
  fi

# install dependencies

  pip --disable-pip-version-check install scipy h5py mhcflurry
  mhcflurry-downloads fetch

  echo "Installing dependencies..."

  cd "$HOME"
  curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz

  Rscript -e \
  'install.packages("devtools"); devtools::install_github("hadley/devtools"); devtools::install_github("r-lib/usethis"); install.packages("testthat")'

  Rscript -e \
  'source("https://bioconductor.org/biocLite.R"); biocLite(c("BiocInstaller", "ShortRead", "biomaRt", "ShortRead"), ask = FALSE)'

  Rscript -e \
  'remove.packages("data.table"); install.packages("data.table", type = "source", repos = "http://Rdatatable.github.io/data.table")'

  Rscript -e \
  'devtools::install_github(c("tidyverse/magrittr", "andrewrech/dt.inflix"))'

# install antigen.garnish

  echo "Installing antigen.garnish..."

##### TODO remove fix_travis ref
  Rscript -e \
  'devtools::install_github("andrewrech/antigen.garnish", ref = "fix_travis")'

# test antigen.garnish

  echo "Testing antigen.garnish installation..."

  mkdir -p "$HOME"/antigen.garnish/tests &&
  cd "$HOME"/antigen.garnish/tests

  Rscript -e \
  'setwd("~/antigen.garnish/tests"); devtools::check("antigen.garnish")'


