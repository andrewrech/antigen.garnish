[![Build Status](http://get.rech.io/build.passing.svg)](http://18.194.224.158:8080/job/antigen.garnish/lastBuild/consoleFull) [![codecov.io](https://codecov.io/github/andrewrech/antigen.garnish/coverage.svg?branch=master)](https://codecov.io/github/andrewrech/antigen.garnish?branch=master) ![](https://img.shields.io/badge/version-0.0.6-blue.svg)

# antigen.garnish

Ensemble neoepitope prediction and multi-parameter quality analysis from direct input, SNVs, indels, and fusions variants in R.

![](http://get.rech.io/antigen.garnish_flowchart.svg)

## Description

An R package for [neoepitope](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, and RNASeq-derived gene fusions and performs ensemble neoepitope prediction using up to 7 algorithms. Input is a VCF file, [JAFFA](https://github.com/Oshlack/JAFFA) output, or table of peptides or transcripts. Outputs are ranked and classified neoepitopes and a summary of neoepitope burden by sample. Neoepitopes are ranked by MHC I/II binding affinity, clonality, RNA expression, dissimilarity to the normal peptidome, and similarity to known immunogenic antigens.

### Advantages

1. **Thoroughness**:
	* missense mutations, insertions, deletions, and gene fusions
	* human and mouse
	* ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
	* ranked by
		* MHC I/II binding affinity
		* clonality
		* RNA expression
		* dissimilarity to the normal peptidome (not presentely in this open source version prior to publication)
		* similarity to known immunogenic antigens
1. **Speed and simplicity**:
	* 1000 variants are ranked in a single step in less than five minutes
	* parallelized using [`parallel::mclapply`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html) and [data.table::setDTthreads](https://github.com/Rdatatable/data.table/wiki), see respective links for information on setting multicore usage
1. **Integration with R/Bioconductor**
	* upstream/VCF processing
	* exploratory data analysis, visualization

## Installation

### Requirements

- Linux
- R &ge; 3.4
- python-pip
- `sudo` is required to install dependencies

or

- a Docker image is available, please contact us

### Install prediction tools and `antigen.garnish`

- detailed installation instructions for bootstrapping a fresh AWS instance can be found in the [wiki](https://github.com/immune-health/antigen.garnish/wiki)

- please note that netMHC, netMHCpan, netMHCII, and netMHCIIpan require academic-use only licenses

One-line complete install command from the shell:

```sh
curl -fsSL http://get.rech.io/antigen.garnish.sh | sudo sh
```

## [Package documentation](https://neoantigens.rech.io/reference/index.html) ([pdf](https://get.rech.io/antigen.garnish.pdf))

### Workflow

  1. Prepare input for MHC affinity prediction and quality analysis:

		* VCF input - `garnish_variants`
		* Fusions from RNASeq via [JAFFA](https://github.com/Oshlack/JAFFA)- `garnish_jaffa`
		* Prepare table of direct transcript or peptide input - see manual page for `?garnish_affinity`

  1. Add MHC alleles of interest - see examples below.
  1. Run ensemble prediction method and perform antigen quality analysis including dissimilarity, IEDB alignment score, and proteome-wide differential agretopicity - `garnish_affinity`.
  1. Summarize output at the sample level with `garnish_summary` and `garnish_plot`, and prioritize the highest quality neoantigens per clone per sample with `garnish_antigens`.

### Examples

#### Predict neoepitopes from missense mutations, insertions, and deletions

```r
library(magrittr)
library(data.table)
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

  # rank results by therapeutic potential
    dt %>%
     garnish_antigens %T>%
      print
```

#### Predict neoepitopes from gene fusions

```r
library(magrittr)
library(data.table)
library(antigen.garnish)

  # load example jaffa output
	dir <- system.file(package = "antigen.garnish") %>%
		file.path(., "extdata/testdata")

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

## Plots and summary tables

- [`garnish_plot`](http://get.rech.io/antigen.garnish_example_plot.pdf)
- [`garnish_antigens`](http://get.rech.io/antigen.garnish_summary_table.png)

## Citation

_Under review._

## Contributing

We welcome contributions and feedback via [Github](https://github.com/immune-health/antigen.garnish/issues) or [email](mailto:leepr@upenn.edu).

## Acknowledgments

We thank the follow individuals for contributions and helpful discussion:

- [David Balli](https://www.linkedin.com/in/davidballi1)
- [Adham Bear](https://www.med.upenn.edu/apps/faculty/index.php/g20001100/p1073)
- [Katelyn Byrne](https://www.parkerici.org/person/katelyn-byrne/)
- [Danny Wells](http://dannykwells.com/)

## License

GNU Lesser General Public License v3.0
