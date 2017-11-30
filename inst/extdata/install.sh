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

  DEBIAN_FRONTEND=noninteractive apt-get -y install ncbi-blast+

  DEBIAN_FRONTEND=noninteractive apt-get -y install gnutls-bin

  DEBIAN_FRONTEND=noninteractive apt-get -y install perlbrew

  cd "/usr/local/bin"
  
  curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz

  wget ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz

  gzip -d Mus_musculus.GRCm38.pep.all.fa.gz

  wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

  gzip -d Homo_sapiens.GRCh38.pep.all.fa.gz

  wget "http://get.rech.io/iedb.fasta"

  makeblastdb -in Mus_musculus.GRCm38.pep.all.fa -input_type fasta -dbtype prot -out mouse.bdb
  makeblastdb -in Homo_sapiens.GRCh38.pep.all.fa -input_type fasta -dbtype prot -out human.bdb
  makeblastdb -in iedb.fasta -input_type fasta -dbtype prot -out iedb.bdb

  Rscript -e \
  'install.packages("devtools", repos = "http://cran.us.r-project.org"); devtools::install_github("hadley/devtools"); install.packages("testthat", repos = "http://cran.us.r-project.org")'

  Rscript -e \
  'source("https://bioconductor.org/biocLite.R"); biocLite(c("ShortRead", "biomaRt", "Biostrings"), suppressUpdates = TRUE, suppressAutoUpdate = TRUE, build_vignettes = FALSE)'

  Rscript -e \
  'install.packages("data.table", type = "source", repos = "http://Rdatatable.github.io/data.table")'

  Rscript -e \
  'devtools::install_github(c("tidyverse/magrittr", "andrewrech/dt.inflix"))'


# install antigen.garnish

  echo "Installing antigen.garnish..."

  Rscript -e \
<<<<<<< HEAD
  'devtools::install_github("andrewrech/antigen.garnish")'
=======
  'devtools::install_github("andrewrech/antigen.garnish@addblast")'
>>>>>>> iedb

  Rscript -e \
  'antigen.garnish::check_pred_tools(); message("Testing antigen.garnish..."); testthat::test_package("antigen.garnish")'
