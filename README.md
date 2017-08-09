[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/antigen.garnish)](http://cran.r-project.org/package=antigen.garnish) ![](https://img.shields.io/badge/build-passing-brightgreen.svg)

# antigen.garnish

Ensemble neoepitope prediction in R.

## Description

An R package for neoepitope analysis that takes human or murine DNA point mutations, insertions, and deletions and performs neoepitope prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/). Output is individual peptides and a summary of neoepitope metrics.

## Installation


### Required prediction tools

Install [mhcflurry](https://github.com/hammerlab/mhcflurry):

```sh
pip install mhcflurry
mhcflurry-downloads fetch
nosetests .
```

Agree to license for non-commercial use and install [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/), and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/). These tools need to be available in `$PATH`.

### R dependencies

Install R dependencies and `antigen.garnish`:

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

* `garnish_variants`: Intake variants from SnpEff.
* `garnish_predictions`: Performs ensemble neoepitope prediction.
* `garnish_summary`: Summarize neoepitope prediction.

### Example and test

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

## Bugs

## Authors

[Andrew J. Rech](mailto:andrewrech@gmail.com)

## License

GNU General Public License v3.0