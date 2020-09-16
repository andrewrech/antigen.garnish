[![rech.io](https://s3.amazonaws.com/get.rech.io/antigen.garnish_build_status.svg)](https://s3.amazonaws.com/get.rech.io/antigen.garnish.test.txt) | [![rech.io](https://img.shields.io/badge/endpoint.svg?url=https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.json)](https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.html) | ![](https://img.shields.io/badge/version-2.0.0-blue.svg) | ![](https://img.shields.io/docker/pulls/leeprichman/antigen_garnish.svg)

# antigen.garnish 2.0

Humand and mouse ensemble tumor neoantigen prediction from SNVs and complex variants. Immunogenicity filtering based on the [Tumor Neoantigen Selection Alliance (TESLA)](https://www.parkerici.org/research-project/tumor-neoantigen-selection-alliance-tesla/).

![](https://get.rech.io/antigen.garnish_flowchart.svg)

## Citation

> Richman LP, Vonderheide RH, and Rech AJ. "Neoantigen dissimilarity to the self-proteome predicts immunogenicity and response to immune checkpoint blockade." Cell Systems. 2019. DOI: [10.1016/j.cels.2019.08.009](https://doi.org/10.1016/j.cels.2019.08.009)

## References

> Rech AJ, Balli D, Mantero A, Ishwaran H, Nathanson KL, Stanger BZ, Vonderheide RH. Tumor immunity and survival as a function of alternative neopeptides in human cancer. Clinical Cancer Research, 2018. DOI: [10.1158/2326-6066.CIR-17-0559](https://cancerimmunolres.aacrjournals.org/content/6/3/276)

> Wells DK, van Buuren MM, Dang KK, Hubbard-Lucey VM, Sheehan KCF, Campbell KM, Lamb A, Ward JP, Sidney J, Blazquez AB, Rech AJ, Zaretsky JM, Comin-Anduix B, Ng AHC, Chour W, Yu TV, Rizvi1 H, Chen JM, Manning P, Steiner GM, Doan XC, The TESLA Consortium, Merghoub T, Guinney J, Kolom A, Selinsky C, Ribas A, Hellmann MD, Hacohen N, Sette A, Heath JR, Bhardwaj N, Ramsdell F, Schreiber RD, Schumacher TN, Kvistborg P, Defranoux N. Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction. Cell. 2020. In press.

## Installation

Two methods exist to run `antigen.garnish`:

1. Docker
2. Linux

### Docker

```sh
docker pull leeprichman/antigen_garnish
```

### Linux

#### Requirements

- R &ge; 3.5.0
- python-pip
- tcsh (required for netMHC)
- `sudo` privileges (required for netMHC)
- GNU Parallel

#### Installation script

The following line downloads and runs the initial [installation script](http://get.rech.io/install_antigen.garnish.sh).

```sh
$ curl -fsSL http://get.rech.io/install_antigen.garnish.sh | sudo sh
```

Next, download the netMHC suite of tools for Linux, available under an academic license:

- [netMHC](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHC)
- [netMHCpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
- [netMHCII](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCII)
- [netMHCIIpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan)

After downloading the files above, move the binaries into the `antigen.garnish` data directory, first setting the `NET_MHC_DIR` and `ANTIGEN_GARNISH_DIR` environmental variables, as shown here:

```sh

NET_MHC_DIR=/path/to/folder/containing/netMHC/downloads
ANTIGEN_GARNISH_DIR=/path/to/antigen.garnish/data/directory

cd "$NET_MHC_DIR" || return 1

mkdir -p "$ANTIGEN_GARNISH_DIR/netMHC" || return 1

find . -name "netMHC*.tar.gz" -exec tar xvzf {} -C "$ANTIGEN_GARNISH_DIR/netMHC" \;

chown "$USER" "$ANTIGEN_GARNISH_DIR/netMHC"
chmod 700 -R "$ANTIGEN_GARNISH_DIR/netMHC"

```

### Development version from master

Follow instructions above under _Installation script_ to install dependencies, and then:

```r
devtools::install_github("immune-health/antigen.garnish")
```

## Package documentation

[Website](https://neoantigens.rech.io/reference/index.html), [PDF](https://get.rech.io/antigen.garnish.pdf).

## Workflow overview

1. Prepare input for MHC affinity prediction and quality analysis:
   - VCF input: (see `?garnish_variants`), or
   - transcript or peptide input (see `?garnish_affinity`)
1. Add MHC alleles of interest (see examples below).
1. Run prediction method (see `?garnish_affinity`)
1. Filter output (see `?garnish_antigens`).

### Running antigen.garnish via Docker

```sh
docker pull leeprichman/antigen_garnish

cID=$(docker run -it -d leeprichman/antigen_garnish /bin/bash)
```

#### Install netMHC binaries (academic license)

Copy netMHC tars to container, changing version names as appropriate:

```sh
netMHC="netMHC-VERSION.Linux.tar.gz"
netMHCII="netMHCII-VERSION.Linux.tar.gz"
netMHCpan="netMHCIIpan-VERSION.Linux.tar.gz"
netMHCIIpan="netMHCpan-VERSION.Linux.tar.gz"

docker cp $netMHC $cID:/$netMHC
docker cp $netMHCII $cID:/$netMHCII
docker cp $netMHCpan $cID:/$netMHCpan
docker cp $netMHCIIpan $cID:/$netMHCIIpan

docker exec $cID config_netMHC.sh
```

#### Interactive use

```sh
docker exec -it $cID bash
```

#### VCF input

Input VCF file name must match the VCF tumor sample name for tumor allelic fraction to be recognized. For many VCFs, the tumor sample name is "TUMOR". In this case, the following input VCF file names would all work:

```
TUMOR.vcf
TUMOR_se.vcf
TUMOR.ann.vcf
TUMOR.vcf.gz
```

MHC input is a JSON file from xHLA or a 2-column tab or comma-separated file ending in "mhc.txt" with the following format:

```
Mouse:
       sample_id,MHC
       mysample.vcf,H-2-Kb H-2-Db H-2-IAb
Human:
       sample_id,MHC
       mysample.vcf,HLA-A*02:01 HLA-B*07:02 or H-2-Kb H-2-Db
```

```sh
VCFO="TUMOR.vcf"
MHC="hla.json"

docker cp $VCFO $cID:/$VCFO
docker cp $MHC $cID:/$MHC

# run antigen.garnish
docker exec $cID run_antigen.garnish.R
```

#### Peptide or transcript input

Start and configure the container as described above. One or more tab or comma-separated input files with the pattern "docker_input" in the file name with the following format are required:

```
Example transcript input

     Column name                 Example input

     sample_id                   sample_1
     ensembl_transcript_id       ENST00000311936
     cDNA_change                 c.718T>A
     MHC                         HLA-A*02:01 HLA-A*03:01


 Example peptide input

     Column name                 Example input

     sample_id                   sample_1
     pep_mut                     MTEYKLVVVDAGGVGKSALTIQLIQNHFV
     mutant_index                25
     MHC                         HLA-A*02:01 HLA-A*03:01
```

```sh
INPUT="docker_input.csv"

docker cp $INPUT $cID:/$DT

# run antigen.garnish
docker exec $cID run_antigen.garnish_direct.R
```

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
    garnish_variants(.) %>%
  # add space separated MHC types
  # see list_mhc() for nomenclature of supported alleles
	# MHC may also be set to "all_human" or "all_mouse" to use all supported alleles
      .[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")] %>%
  # predict neoantigens
    garnish_affinity(.) %>%
```

#### Directly calculate IEDB score and dissimilarity for a list of sequences

```r
library(magrittr)
library(data.table)
library(antigen.garnish)

  # generate our character vector of sequences
  v <- c("SIINFEKL", "ILAKFLHWL", "GILGFVFTL")
  # calculate IEDB score
  v %>% iedb_score(db = "human") %>% print
	# calculate dissimilarity
	v %>% garnish_dissimilarity(db = "human") %>% print
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

## Contributing

We welcome contributions and feedback via [Github](https://github.com/immune-health/antigen.garnish/issues) or [email](mailto:leepr@upenn.edu).

## Acknowledgments

We thank the follow individuals for contributions and helpful discussion:

- [David Balli](https://www.linkedin.com/in/davidballi1)
- [Adham Bear](https://www.med.upenn.edu/apps/faculty/index.php/g20001100/p1073)
- [Katelyn Byrne](https://www.parkerici.org/person/katelyn-byrne/)
- [Danny Wells](http://dannykwells.com/)

## License

Please see [LICENSE](https://github.com/immune-health/antigen.garnish/blob/master/LICENSE).
