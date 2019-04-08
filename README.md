[![Build Status](http://get.rech.io/build.passing.svg)](http://18.194.224.158:8080/job/antigen.garnish/lastBuild/consoleFull) [![codecov.io](https://codecov.io/github/andrewrech/antigen.garnish/coverage.svg?branch=master)](https://codecov.io/github/andrewrech/antigen.garnish?branch=master) ![](https://img.shields.io/badge/version-0.0.6-blue.svg)

# antigen.garnish

Ensemble neoepitope prediction from DNA variants in R.

![](http://get.rech.io/antigen.garnish_flowchart.svg)

## Description

An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, and gene fusions and performs ensemble neoepitope prediction using 7 algorithms. Input is a VCF file or table of peptides. Output is ranked neoepitopes and a summary of neoepitope load and fitness by sample. Neoepitopes are ranked by MHC I/II binding affinity, clonality, RNA expression, dissimilarity to the normal peptidome, and similarity to known immunogenic antigens. [Package schematic](https://github.com/andrewrech/antigen.garnish/wiki/Package-schematic). [More information.](http://antigen-garnish-presentation.s3-website-us-east-1.amazonaws.com)

### Advantages

1. **Thoroughness**:
		* missense mutations, insertions, deletions, and gene fusions
		* human and mouse
		* ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
		* ranked by
			* MHC I/II binding affinity
			* clonality
			* RNA expression
			* dissimilarity to the normal peptidome
			* similarity to known immunogenic antigens
1. **Speed and simplicity**:
		* 1000 variants are ranked in a single step in less than five minutes
1. **Integration with R/Bioconductor**
		* upstream/VCF processing
		* exploratory data analysis, visualization

## Installation

### Requirements

- Linux
- R &ge; 3.4
- python-pip

### Install prediction tools and `antigen.garnish`

```sh
curl -fsSL http://get.rech.io/antigen.garnish.sh | sudo sh
```

- detailed installation instructions can be found in the [wiki](https://github.com/andrewrech/antigen.garnish/wiki).

## Package functions

- `garnish_variants`: process missense / indel VCF variants from [SnpEff](http://snpeff.sourceforge.net/)
- `garnish_jaffa`: process gene fusions from [JAFFA](https://github.com/Oshlack/JAFFA)
- `garnish_affinity`: perform ensemble neoepitope prediction
- `garnish_summary`: summarize and rank results at sample-level
- `garnish_antigens`: list top neoepitope sequences, sources, and properties for each[clone](https://github.com/lima1/PureCN) and sample
- `garnish_plot`: generate sample-level summary plots
- `garnish_slim`: reduce an output table to fewest columns containing minimal prediction information

### Examples

#### Predict neoepitopes from missense mutations, insertions, and deletions

```r
library(magrittr)
library(antigen.garnish)

  # load an example VCF
		dir <- system.file(package = "antigen.garnish") %>%
			file.path(., "extdata/testdata")

    dt <- "antigen.garnish_example.vcf" %>%
    file.path(dir, .) %>%

  # extract variants
    garnish_variants %>%

  # add space separated MHC types
  # see list_mhc() for nomenclature of supported alleles
  # separate murine and human alleles into separate rows, even if same sample_id.

      .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
                   "H-2-Kb H-2-IAd",
                   "HLA-A*01:47 HLA-DRB1*03:08")] %>%

  # predict neoepitopes
    garnish_affinity

  # summarize predictions
    dt %>%
      garnish_summary %T>%
        print

  # generate summary graphs
    dt %>% garnish_plot
```

#### Predict neoepitopes from gene fusions

```r
library(magrittr)
library(antigen.garnish)

	dir <- system.file(package = "antigen.garnish") %>%
		file.path(., "extdata/testdata")

  # load example jaffa output
  path <- "antigen.garnish_jaffa_results.csv" %>%
      file.path(dir, .)
  fasta_path <- "antigen.garnish_jaffa_results.fasta" %>%
      file.path(dir, .)

  # get predictions
    dt <- garnish_jaffa(path, db = "GRCm38", fasta_path) %>%

  # add MHC info with list_mhc() compatible names
    .[, MHC := "H-2-Kb"] %>%

  # get predictions
    garnish_affinity %>%

  # summarize predictions
    garnish_summary %T>%
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
    make_nmers %T>% print
```

## Bugs

## Authors

- [Lee P. Richman](http://www.med.upenn.edu/apps/faculty/index.php/g275/p1073)
- [David Balli](https://www.linkedin.com/in/davidballi1)
- [Robert H. Vonderheide](https://www.med.upenn.edu/apps/faculty/index.php/g20000320/p1073)
- [Andrew J. Rech](https://rech.io)

## Citation

Rech AJ, Balli D, Stanger BZ, Vonderheide RH. Tumor immunity and survival as a function of alternative neoepitopes in human cancer. Cancer Immunology Research, 2018 Jan 16. pii: canimm.0559.2017. PMID: [29339376](https://www.ncbi.nlm.nih.gov/pubmed/29339376)

## Contributing

We welcome contributions and feedback via Github or [email](mailto:rech@rech.io).

## License

GNU General Public License v3.0
