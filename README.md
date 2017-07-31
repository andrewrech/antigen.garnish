[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/antigen.garnish)](http://cran.r-project.org/package=antigen.garnish) ![](https://img.shields.io/badge/build-passing-brightgreen.svg)

# antigen.garnish

Ensemble neoepitope prediction in R.

## Description

An R package for neoepitope analysis that takes human or murine missense variants and returns a summary of ensemble neoepitope prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/).

## Installation

```r
devtools::install_github("andrewrech/dt.inflix")
devtools::install_github("andrewrech/antigen.garnish")
```

### Requirements

* antigen prediction tools in `$PATH`

## Manifest

* `garnish_variants`: Intake variants from SnpEff.
* `garnish_predictions`: Performs ensemble epitope prediction.
* `garnish_summary`: Summarize epitope prediction.

### Example

```r
library(testthat)
library(antigen.garnish)
library(data.table)

mhc_dt <- data.table::data.table(
            transcript_affected = c("ENST00000256078", "ENST00000256078"),
            sample_id = c("test_sample_1", "test_sample_2"),
            aa_mutation = c("G12D", "G13D"),
            MHC = c("HLA-A*02:01 HLA-A*03:01 HLA-DRB10301 HLA-DRB1*14:67",
                    "HLA-A*02:01 HLA-A*03:01 HLA-DRB1*03:01 HLA-DRB1*14:67"))

mhc_dt <- garnish_predictions(mhc_dt)

dt <- garnish_summary(mhc_dt)

testthat::all.equal(dt,
    structure(list(
     sample_id = c("test_sample_1", "test_sample_2"),
     ag_cdn = c(0L, 0L),
     ag_adn = c(0L, 0L),
     ag_adn_top = c(5.20144265891428, 1.56599284273263),
     ag_cdn_top = c(0.00954700205453025, 0.00954700205453025),
     ag_binders = c(11L, 9L),
     ag_peptides = c(74L, 76L),
     ag_predictions = c(148L, 152L)),
       .Names = c("sample_id", "ag_cdn", "ag_adn", "ag_adn_top",
       "ag_cdn_top", "ag_binders", "ag_peptides", "ag_predictions"),
       row.names = c(NA, -2L),
       class = c("data.table", "data.frame"))
          )

```

## Bugs

## Authors

[Andrew J. Rech](mailto:andrewrech@gmail.com)

## License

GNU General Public License v3.0