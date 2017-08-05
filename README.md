[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/antigen.garnish)](http://cran.r-project.org/package=antigen.garnish) ![](https://img.shields.io/badge/build-passing-brightgreen.svg)

# antigen.garnish

Ensemble neoepitope prediction in R.

## Description

An R package for neoepitope analysis that takes human or murine missense variants and returns a summary of ensemble neoepitope prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/).

## Installation

Install [mhcflurry](https://github.com/hammerlab/mhcflurry):

```sh
pip install mhcflurry
mhcflurry-downloads fetch
nosetests .
```

Agree to license for non-commercial use and install [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/), and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/). These tools need to be available in `$PATH`.

Install R dependencies and `antigen.garnish`:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

devtools::install_github("tidyverse/magrittr")
devtools::install_github("andrewrech/dt.inflix")

devtools::install_github("andrewrech/antigen.garnish")
```

## [Package documentation](http://get.rech.io/antigen.garnish.pdf)

* `garnish_variants`: Intake variants from SnpEff.
* `garnish_predictions`: Performs ensemble neoepitope prediction.
* `garnish_summary`: Summarize neoepitope prediction.

```r
# Generate package documentation

system(paste(shQuote(file.path(R.home("bin"), "R")),
    "CMD", "Rd2pdf", shQuote(find.package("antigen.garnish"))))

```

### Example

```r
library(magrittr)

    # download an example VCF
    "antigen.garnish_example.vcf" %T>%
    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%

    # extract variants
    antigen.garnish::garnish_variants %>%

    # add MHC types
    .$antigen.garnish_input %>%
        .[, MHC := c("HLA-A*02:01 HLA-A*03:01 HLA-DRB10301 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67")] %>%

    # predict
    antigen.garnish::garnish_predictions %>%

    # summarize
    antigen.garnish::garnish_summary %T>%

    print %>%

    # does antigen.garnish work?
    testthat::compare(
            data.table::fread("http://get.rech.io/ag_test.csv"))
```

## Bugs

## Authors

[Andrew J. Rech](mailto:andrewrech@gmail.com)

## License

GNU General Public License v3.0