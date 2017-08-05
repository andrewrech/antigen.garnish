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


```r
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

devtools::install_github("tidyverse/magrittr")
devtools::install_github("andrewrech/dt.inflix")

devtools::install_github("andrewrech/antigen.garnish")
```

## Manifest

* `garnish_variants`: Intake variants from SnpEff.
* `garnish_predictions`: Performs ensemble neoepitope prediction.
* `garnish_summary`: Summarize neoepitope prediction.

### Example

```r
library(magrittr)

dt <-
    # load an example VCF
    system.file("extdata",
          "antigen.garnish_example.vcf",
          package = "antigen.garnish") %>%

    # extract variants
    antigen.garnish::garnish_variants(.)


    # add MHC types
    package_test <- dt$antigen.garnish_input %>%
        .[, MHC := c("HLA-A*02:01 HLA-A*03:01 HLA-DRB10301 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67")] %>%

    # predict
    antigen.garnish::garnish_predictions(.) %>%

    # sumarize
    antigen.garnish::garnish_summary(.)


    # does antigen.garnish work?
    testthat::compare(package_test,
    structure(list(sample_id = "tumor",
                    priority_neos = 0L,
                    classic_neos = 0L,
                    alt_neos = 2L,
                    alt_neos_top = 32.5590961308976,
                    classic_neos_top = 0.0086649901329808,
                    binders = 7L,
                    peptides = 231L,
                    predictions = 462L),
                    .Names = c("sample_id",
                    "priority_neos", "classic_neos", "alt_neos", "alt_neos_top",
                    "classic_neos_top", "binders", "peptides", "predictions"),
                    row.names = c(NA, -1L),
                    class = c("data.table", "data.frame"))
    )
```

## Bugs

## Authors

[Andrew J. Rech](mailto:andrewrech@gmail.com)

## License

GNU General Public License v3.0