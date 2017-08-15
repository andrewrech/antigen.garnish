[![Build Status](https://travis-ci.org/andrewrech/antigen.garnish.svg?branch=master)](https://travis-ci.org/andrewrech/antigen.garnish) [![codecov.io](https://codecov.io/github/andrewrech/antigen.garnish/coverage.svg?branch=master)](https://codecov.io/github/andrewrech/antigen.garnish?branch=master) ![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/antigen.garnish) ![](https://img.shields.io/badge/version-0.0.1-blue.svg)



# antigen.garnish

Ensemble neoepitope prediction from DNA variants in R.

![](http://get.rech.io/antigen.garnish_flowchart.svg)

## Description

An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA point mutations, insertions, and deletions in VCF format and performs neoepitope prediction. Output is individual peptides and a summary by sample.

### Advantages

1. **Simplicity**: summarized neoepitopes for each sample
1. **Thoroughness**:
    - missense mutations and frameshifts
    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry),
    [mhcnuggets](https://github.com/KarchinLab/mhcnuggets),
    [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/), and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i).
1. **Speed**:
    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
    - on an Amazon Web Services `m4.16xlarge` EC2 instance, 20,000 consensus predictions using 100+ MHC types in under 5 minutes

## Installation

### Requirements

* Linux
* R &ge; 3.4

### Install required prediction tools

Install [mhcflurry](https://github.com/hammerlab/mhcflurry):

```sh
pip install mhcflurry
mhcflurry-downloads fetch
nosetests .
```

Install [mhcnuggets](https://github.com/KarchinLab/mhcnuggets):

```sh
cd "$HOME"
git clone https://github.com/KarchinLab/mhcnuggets.git
```

Install [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/), and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/):

```sh
  cd "$HOME"
  wget "http://get.rech.io/netMHC.tar.gz"
  tar -xvzf netMHC.tar.gz
```

### Install antigen.garnish

```r
devtools::install_github("andrewrech/antigen.garnish")
```

## Package documentation

* `garnish_variants`: process variants from [SnpEff](http://snpeff.sourceforge.net/)
* `garnish_predictions`: perform ensemble neoepitope prediction
* `garnish_summary`: summarize neoepitope prediction

### Generate documentation

```r
system(paste(shQuote(file.path(R.home("bin"), "R")),
    "CMD", "Rd2pdf", shQuote(find.package("antigen.garnish"))))
```

### Examples

#### Predict neoepitopes

```r
library(magrittr)

    # download an example VCF
    dt <- "antigen.garnish_example.vcf" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

    # extract variants
    antigen.garnish::garnish_variants %>%

    # add MHC types
        .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-DRB1*14:67",
                    "HLA-A*03:01 HLA-DRB1*03:01")] %>%

    # predict neoepitopes
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
       mutant_loc = c(5, 25, 47, 50),
       pep_type = "test",
       var_uuid = c("middle",
                    "back_truncate",
                    "front_truncate",
                    "end")) %>%
    # create nmers
    antigen.garnish::garnish_predictions_worker %T>% print
```

## Bugs

## Authors

* [Andrew J. Rech](http://info.rech.io)
* [Robert H. Vonderheide](http://www.med.upenn.edu/apps/faculty/index.php/g275/p1073)

## License

GNU General Public License v3.0
