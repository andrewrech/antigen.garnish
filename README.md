[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/antigen.garnish)](http://cran.r-project.org/package=antigen.garnish)

# antigen.garnish

Ensemble neoepitope prediction from DNA variants in R.

![](http://get.rech.io/antigen.garnish_flowchart.svg)

## Description

An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA point mutations, insertions, and deletions in VCF format and performs neoepitope prediction. Output is individual peptides and a summary by sample.

### Advantages

1. **Simplicity**: summarized neoepitopes for each sample
1. **Thoroughness**:
    - missense mutations and frameshifts
    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/).
1. **Speed**:
    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
    - parallelized prediction

## Installation


### Required prediction tools

Install [mhcflurry](https://github.com/hammerlab/mhcflurry):

```sh
pip install mhcflurry
mhcflurry-downloads fetch
nosetests .
```

Install [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/), and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/). These tools need to be available in `$PATH`.

### R dependencies

Install R dependencies:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

devtools::install_github("tidyverse/magrittr")  ## Github version is required
devtools::install_github("andrewrech/dt.inflix")
```

### antigen.garnish

```r
devtools::install_github("andrewrech/antigen.garnish")
```

## [Package documentation](http://get.rech.io/antigen.garnish.pdf)

* `garnish_variants`: Process variants from [SnpEff](http://snpeff.sourceforge.net/).
* `garnish_predictions`: Perform ensemble neoepitope prediction.
* `garnish_summary`: Summarize neoepitope prediction.

### Examples

#### Predict neoepitopes

```r
library(magrittr)

    # download an example VCF
    "antigen.garnish_example.vcf" %T>%
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

    print %>%

    # does antigen.garnish work?
    testthat::compare(
       data.table::fread("http://get.rech.io/antigen.garnish_example_summary.txt"))
```

#### Test the sliding window

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