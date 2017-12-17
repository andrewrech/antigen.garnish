[![Build Status](https://travis-ci.org/andrewrech/antigen.garnish.svg?branch=master)](https://travis-ci.org/andrewrech/antigen.garnish) [![codecov.io](https://codecov.io/github/andrewrech/antigen.garnish/coverage.svg?branch=master)](https://codecov.io/github/andrewrech/antigen.garnish?branch=master) ![](https://img.shields.io/badge/version-0.0.5-blue.svg)



# antigen.garnish

Ensemble neoepitope prediction from DNA variants in R.

![](http://get.rech.io/antigen.garnish_flowchart.svg)

## Description

An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, and gene fusions and performs neoepitope prediction. Input is a VCF file or table of peptides. Output is neoepitopes and a summary of neoepitope load and fitness by sample. [More information.](http://antigen-garnish-presentation.s3-website-us-east-1.amazonaws.com)

### Advantages

1. **Simplicity**:
    - VCF or table input
    - summarized neoepitopes for each sample
1. **Thoroughness**:
    - missense mutations, insertions, deletions, and gene fusions
    - ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
    - filter against all known normal proteins and rank by fitness
1. **Speed**:
    - produce all possible 8-15-mer peptides from 10,000 variants in under 1 minute on a normal laptop
    - on an Amazon Web Services `m4.16xlarge` EC2 instance, 20,000 consensus predictions using 100+ MHC types in under 5 minutes

## Installation

### Requirements

* Linux
* R &ge; 3.4
* python-pip

### Install prediction tools and `antigen.garnish`

```sh
curl -fsSL http://get.rech.io/install_iedb_branch.sh | sudo sh
```

## [Package documentation](http://get.rech.io/antigen.garnish.pdf)

* `garnish_variants`: process VCF variants from [SnpEff](http://snpeff.sourceforge.net/)
* `garnish_jaffa`: process gene fusions from [JAFFA](https://github.com/Oshlack/JAFFA)
* `garnish_predictions`: perform ensemble neoepitope prediction
* `garnish_summary`: summarize neoepitope prediction
* `garnish_plot`: generate summary plots for neoepitope prediction
* `list_mhc`: list all supported MHC allele syntax

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

  # add space separated MHC types
  # see antigen.garnish::list_mhc() for nomenclature of supported alleles

      .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                   "H-2-Kb H-2-IAd",
                   "HLA-A*01:47 HLA-DRB1*03:08")] %>%

  # predict neoepitopes
    antigen.garnish::garnish_predictions

  # summarize predictions
    dt %>%
      antigen.garnish::garnish_summary %T>%
        print

  # generate summary graphs
    dt %>% garnish_plot

  # apply fitness model from Luksza et al.
    dt %>% garnish_fitness

```

#### Predict neoepitopes from gene fusions

```r
library(magrittr)
library(antigen.garnish)

  # load example jaffa output
    path <- "antigen.garnish_jaffa_results.csv" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.csv", .)
    fasta_path <- "antigen.garnish_jaffa_results.fasta" %T>%
      utils::download.file("http://get.rech.io/antigen.garnish_jaffa_results.fasta", .)

  # get predictions
    dt <- antigen.garnish::garnish_jaffa(path, db = "GRCm38", fasta_path) %>%

  # add MHC info with antigen.garnish::list_mhc() compatible names
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

* [Lee P. Richman](http://www.med.upenn.edu/apps/faculty/index.php/g275/p1073)
* [David Balli](https://www.linkedin.com/in/davidballi1)
* [Robert H. Vonderheide](https://www.med.upenn.edu/apps/faculty/index.php/g20000320/p1073)
* [Andrew J. Rech](http://rech.io)

## Contributing

We welcome contributions and feedback via Github or [email](mailto:andrewrech@gmail.com).

## License

GNU General Public License v3.0
