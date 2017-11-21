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

  echo "Installing dependencies..."

  pip --disable-pip-version-check install scipy h5py mhcflurry biopython
  mhcflurry-downloads fetch

  sudo DEBIAN_FRONTEND=noninteractive apt-get -y install ncbi-blast+

  sudo DEBIAN_FRONTEND=noninteractive apt-get -y install gnutls-bin

  cd "/usr/bin"
  sudo curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | sudo tar -xvz

  Rscript -e \
  'install.packages("devtools", repos="http://cran.us.r-project.org"); devtools::install_github("hadley/devtools"); install.packages("testthat", repos="http://cran.us.r-project.org")'

  Rscript -e \
  'source("https://bioconductor.org/biocLite.R"); biocLite(c("ShortRead", "biomaRt", "Biostrings"), ask = FALSE)'

  Rscript -e \
  'install.packages("data.table", type = "source", repos = "http://Rdatatable.github.io/data.table")'

  Rscript -e \
  'devtools::install_github(c("tidyverse/magrittr", "andrewrech/dt.inflix"))'


# install antigen.garnish

  echo "Installing antigen.garnish..."

  Rscript -e \
  'devtools::install_github("andrewrech/antigen.garnish@addblast")'

  Rscript -e \
  'antigen.garnish::check_pred_tools()'

  echo "Testing antigen.garnish..."

  Rscript -e \
  'testthat::test_package("antigen.garnish")'
