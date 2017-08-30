[![Build Status](https://travis-ci.org/andrewrech/antigen.garnish.svg?branch=master)](https://travis-ci.org/andrewrech/antigen.garnish) [![codecov.io](https://codecov.io/github/andrewrech/antigen.garnish/coverage.svg?branch=master)](https://codecov.io/github/andrewrech/antigen.garnish?branch=master) ![](https://img.shields.io/badge/version-0.0.1-blue.svg)



# antigen.garnish

Ensemble neoepitope prediction from DNA variants in R.

![](http://get.rech.io/antigen.garnish_flowchart.svg)

## Description

An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA point mutations, insertions, and deletions in VCF format and performs neoepitope prediction. Output is individual peptides and a summary by sample.

### Advantages

1. **Simplicity**:
    - table or vcf input
    - summarized neoepitopes for each sample
1. **Thoroughness**:
    - missense mutations and frameshifts
    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i).
1. **Speed**:
    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
    - on an Amazon Web Services `m4.16xlarge` EC2 instance, 20,000 consensus predictions using 100+ MHC types in under 5 minutes

## Installation

### Requirements

* Linux
* R &ge; 3.4

### Install required prediction tools

Install `mhcflurry` and `mhcnuggets` dependencies from the command line:

```sh
pip install mhcflurry scipy h5py
mhcflurry-downloads fetch
nosetests .
```

Install `netMHC` and `mhcnuggets` prediction tools from the command line:

```sh
  cd "$HOME"
  wget "http://get.rech.io/antigen.garnish.tar.gz"
  tar -xvzf antigen.garnish.tar.gz
```

### Install antigen.garnish

From R:

```r
if (!"devtools" %in% installed.packages()) install.packages("devtools")

devtools::install_github("andrewrech/antigen.garnish")
```

## Package [documentation](http://get.rech.io/antigen.garnish.pdf)

* `garnish_variants`: process variants from [SnpEff](http://snpeff.sourceforge.net/)
* `garnish_predictions`: perform ensemble neoepitope prediction
* `garnish_summary`: summarize neoepitope prediction

### Examples

#### Predict neoepitopes

```r
library(magrittr)

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
       var_uuid = c("middle",
                    "back_truncate",
                    "front_truncate",
                    "end")) %>%
    # create nmers
    antigen.garnish::get_nmers %T>% print
```

## Bugs

## Authors

* [Andrew J. Rech](http://info.rech.io)
* [Lee P. Richman](http://www.med.upenn.edu/apps/faculty/index.php/g275/p1073)
* [Robert H. Vonderheide](http://www.med.upenn.edu/apps/faculty/index.php/g275/p1073)

## License

GNU General Public License v3.0
