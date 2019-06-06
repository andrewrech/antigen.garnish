[![rech.io](https://s3.amazonaws.com/get.rech.io/antigen.garnish_build_status.svg)](https://s3.amazonaws.com/get.rech.io/antigen.garnish.test.txt) | [![rech.io](https://img.shields.io/badge/endpoint.svg?url=https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.json)](https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.html) | ![](https://img.shields.io/badge/version-1.0.0-blue.svg) | ![](https://img.shields.io/docker/pulls/leeprichman/antigen_garnish.svg)

# antigen.garnish

Ensemble tumor neoantigen prediction and multi-parameter quality analysis from direct input, SNVs, indels, or gene fusion variants.

![](https://get.rech.io/antigen.garnish_flowchart.svg)

[Detailed flowchart.](https://get.rech.io/antigen.garnish_flowchart_detailed.svg)

## Description

An R package for [neoantigen](http://science.sciencemag.org/content/348/6230/69) analysis that takes human or murine DNA missense mutations, insertions, deletions, or RNASeq-derived gene fusions and performs ensemble neoantigen prediction using 7 algorithms. Input is a VCF file, [JAFFA](https://github.com/Oshlack/JAFFA) output, or table of peptides or transcripts. Outputs are ranked and summarized by sample. Neoantigens are ranked by MHC I/II binding affinity, clonality, RNA expression, similarity to known immunogenic antigens, and dissimilarity to the normal peptidome.

### Advantages

1. **Thoroughness**:
	* missense mutations, insertions, deletions, and gene fusions
	* human and mouse
	* ensemble MHC class I/II binding prediction using [mhcflurry](https://github.com/hammerlab/mhcflurry), [mhcnuggets](https://github.com/KarchinLab/mhcnuggets-2.0), [netMHC](http://www.cbs.dtu.dk/services/NetMHC/), [netMHCII](http://www.cbs.dtu.dk/services/NetMHCII/), [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) and [netMHCIIpan](http://www.cbs.dtu.dk/services/NetMHCIIpan/i)
	* ranked by
		* MHC I/II binding affinity
		* clonality
		* RNA expression
		* similarity to known immunogenic antigens
		* dissimilarity to the normal peptidome
2. **Speed and simplicity**:
	* 1000 variants are ranked in a single step in less than five minutes
	* parallelized using [`parallel::mclapply`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html) and [data.table::setDTthreads](https://github.com/Rdatatable/data.table/wiki), see respective links for information on setting multicore usage
3. **Integration with R/Bioconductor**
	* upstream/VCF processing
	* exploratory data analysis, visualization

## Installation

### Requirements

- Linux
- R &ge; 3.4
  - see [documentation](https://get.rech.io/antigen.garnish.pdf) `Imports` for R, [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) package from Bioconductor
- python-pip
- `sudo` is required to install prediction tool dependencies

or

- a Docker image is available, please see the [wiki](https://github.com/immune-health/antigen.garnish/wiki) or [contact us](mailto:leepr@upenn.edu)

### Install all dependencies, prediction tools, and `antigen.garnish`

One-line [installation script](http://get.rech.io/install_antigen.garnish.sh):

```sh
$ curl -fsSL http://get.rech.io/install_antigen.garnish.sh | sudo sh
```

- if installing without using the above installation script, set `$AG_DATA_DIR` to the [required data directory](http://get.rech.io/antigen.garnish.tar.gz):

```sh
$ curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz
$ export AG_DATA_DIR="$PWD/antigen.garnish"
```

- detailed installation instructions for bootstrapping a fresh AWS instance can be found in the [wiki](https://github.com/immune-health/antigen.garnish/wiki)
- please note that netMHC, netMHCpan, netMHCII, and netMHCIIpan require academic-use only licenses

## [Package documentation](https://neoantigens.rech.io/reference/index.html) ([pdf](https://get.rech.io/antigen.garnish.pdf))

### Workflow

  1. Prepare input for MHC affinity prediction and quality analysis:

		* VCF input - `garnish_variants`
		* Fusions from RNASeq via [JAFFA](https://github.com/Oshlack/JAFFA)- `garnish_jaffa`
		* Prepare table of direct transcript or peptide input - see manual page in R (`?garnish_affinity`)

  1. Add MHC alleles of interest - see examples below.
  1. Run ensemble prediction method and perform antigen quality analysis including proteome-wide differential agretopicity, IEDB alignment score, and dissimilarity: `garnish_affinity`.
  1. Summarize output by sample level with `garnish_summary` and `garnish_plot`, and prioritize the highest quality neoantigens per clone and sample with `garnish_antigens`.

### Examples

#### Predict neoantigens from missense mutations, insertions, and deletions

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

      .[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")] %>%

  # predict neoantigens
    garnish_affinity

  # summarize predictions
    dt %>%
      garnish_summary %T>%
        print

  # generate summary graphs
    dt %>% garnish_plot
```

#### Predict neoantigens from gene fusions

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

#### Get full MHC affinity output from a Microsoft Excel file of variants

```r
library(magrittr)
library(data.table)
library(antigen.garnish)

  # load example Microsoft Excel file
  dir <- system.file(package = "antigen.garnish") %>%
    file.path(., "extdata/testdata")

  path <- "antigen.garnish_test_input.xlsx" %>%
    file.path(dir, .)

  # predict neoantigens
    dt <- garnish_affinity(path = path) %T>%
      str
```

#### Automated testing

From ./`<Github repo>`:

```r
  devtools::test(reporter = "summary")
```

#### How are peptides generated?

```r
  library(magrittr)
  library(data.table)
  library(antigen.garnish)

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

* `garnish_plot` output:

![](https://get.rech.io/antigen.garnish_example_plot.png)

* `garnish_antigens` output:

![](https://get.rech.io/antigen.garnish_summary_table.png)

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
