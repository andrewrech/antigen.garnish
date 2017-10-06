[![Build Status](https://travis-ci.org/andrewrech/antigen.garnish.svg?branch=master)](https://travis-ci.org/andrewrech/antigen.garnish) [![codecov.io](https://codecov.io/github/andrewrech/antigen.garnish/coverage.svg?branch=master)](https://codecov.io/github/andrewrech/antigen.garnish?branch=master) ![](https://img.shields.io/badge/version-0.0.4-blue.svg)



# antigen.garnish

Ensemble neoepitope prediction from DNA variants in R.

![](http://get.rech.io/antigen.garnish_flowchart.svg)

## Description

An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, and gene fusions and performs neoepitope prediction. Input is a VCF file or table of peptides. Output is neoepitopes and a summary by sample. [More information.](http://antigen-garnish-presentation.s3-website-us-east-1.amazonaws.com)

### Advantages

1. **Simplicity**:
    - VCF or table input
    - summarized neoepitopes for each sample
1. **Thoroughness**:
    - missense mutations, insertions, deletions, and gene fusions
    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i).
    - filter against all known normal proteins
1. **Speed**:
    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
    - on an Amazon Web Services `m4.16xlarge` EC2 instance, 20,000 consensus predictions using 100+ MHC types in under 5 minutes

## Installation

### Requirements

* Linux
* R &ge; 3.4

### Install depdendencies

Install R (version 3.4), prediction tools, and dependencies on a Ubuntu [AMI](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html) (Ubuntu Server 16.04 LTS (HVM) - ami-cd0f5cb6):

```sh
cd "$HOME"
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 &&
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
sudo apt-get update -y
sudo apt-get install -y r-base python-pip libcurl4-gnutls-dev libssl-dev subversion libxml2-dev

sudo pip install scipy h5py mhcflurry
mhcflurry-downloads fetch

wget "http://get.rech.io/antigen.garnish.tar.gz"
tar -xvzf antigen.garnish.tar.gz
```

Install antigen.garnish and dependencies from R:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
install.packages("devtools")
devtools::install_github("andrewrech/antigen.garnish")
```

## [Package documentation](http://get.rech.io/antigen.garnish.pdf)

* `garnish_variants`: process VCF variants from [SnpEff](http://snpeff.sourceforge.net/)
* `garnish_jaffa`: process gene fusions from [JAFFA](https://github.com/Oshlack/JAFFA)
* `garnish_predictions`: perform ensemble neoepitope prediction
* `garnish_summary`: summarize neoepitope prediction

### Examples

#### Predict neoepitopes from missense mutations, insertions, and deletions

```r
library(magrittr)
library(antigen.garnish)

  # download an example VCF
    dt <- "antigen.garnish_example.vcf" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

  # extract variants
    antigen.garnish::garnish_variants %>%

  # add test MHC types
      .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                   "H-2-Kb H-2-IAd",
                  "HLA-A*01:47 HLA-DRB1*03:08")] %>%

  # predict neoepitopes
    antigen.garnish::garnish_predictions %>%

  # summarize predictions
    antigen.garnish::garnish_summary %T>%
    print
```
#### Predict neoepitopes from gene fusions

```r
library(magrittr)
library(antigen.garnish)

  # load some test jaffa output data
    path <- "antigen.garnish_jaffa_results.csv" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.csv", .)
    fasta_path <- "antigen.garnish_jaffa_results.fasta" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)

  # get predictions
    dt <- antigen.garnish::garnish_jaffa(path, db = "GRCm38", fasta_path) %>%

  # add MHC info
    .[, MHC := "H-2-Kb"] %>%

  # get predictions
    antigen.garnish::garnish_predictions %>%

  # summarize predictions
    antigen.garnish::garnish_summary %T>%
    print
```

### Tests

#### Automated testing

```r
library(testthat)

testthat::test_package("antigen.garnish")
```

#### How are peptides generated?

```r
library(magrittr)

  # generate a fake peptide
    dt <- data.table::data.table(
       pep_base = "Y___*___THIS_IS_________*___A_CODE_TEST!______*__X",
       mutant_index = c(5, 25, 47, 50),
       pep_type = "test",
       var_uuid = c(
                    "front_truncate",
                    "middle",
                    "back_truncate",
                    "end")) %>%
  # create nmers
    antigen.garnish::make_nmers %T>% print
```

## Bugs

## Authors

* [Andrew J. Rech](http://info.rech.io) (maintainer)
* [Lee P. Richman](http://www.med.upenn.edu/apps/faculty/index.php/g275/p1073)
* [David Balli](https://www.linkedin.com/in/davidballi1)
* [Robert H. Vonderheide](http://www.med.upenn.edu/apps/faculty/index.php/g275/p1073)

## Contributing

We welcome contributions and feedback via Github or [email](mailto:andrewrech@gmail.com).

## License

GNU General Public License v3.0
