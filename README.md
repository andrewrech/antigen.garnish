[![rech.io](https://s3.amazonaws.com/get.rech.io/antigen.garnish_build_status.svg)](https://s3.amazonaws.com/get.rech.io/antigen.garnish.test.txt) | [![rech.io](https://img.shields.io/badge/endpoint.svg?url=https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.json)](https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.html) | ![](https://img.shields.io/badge/version-2.0.0-blue.svg) | ![](https://img.shields.io/docker/pulls/leeprichman/antigen_garnish.svg)

# antigen.garnish 2.0

Human and mouse ensemble tumor neoantigen prediction from SNVs and complex variants. Immunogenicity filtering based on the [Tumor Neoantigen Selection Alliance (TESLA)](https://www.parkerici.org/research-project/tumor-neoantigen-selection-alliance-tesla/).

![](https://get.rech.io/antigen.garnish_flowchart.svg)

## Citation

> Richman LP, Vonderheide RH, and Rech AJ. Neoantigen dissimilarity to the self-proteome predicts immunogenicity and response to immune checkpoint blockade. Cell Systems. 2019.

## Selected references

> Duan, F., Duitama, J., Seesi, S.A., Ayres, C.M., Corcelli, S.A., Pawashe, A.P., Blanchard, T., McMahon, D., Sidney, J., Sette, A., et al. Genomic and bioinformatic profiling of mutational neoepitopes reveals new rules to predict anticancer immunogenicity. J Exp Med. 2014.

> Luksza, M, Riaz, N, Makarov, V, Balachandran VP, et al. A neoepitope fitness model predicts tumour response to checkpoint blockade immunotherapy. Nature. 2017.

> Rech AJ, Balli D, Mantero A, Ishwaran H, Nathanson KL, Stanger BZ, Vonderheide RH. Tumor immunity and survival as a function of alternative neopeptides in human cancer. Clinical Cancer Research, 2018.

> Wells DK, van Buuren MM, Dang KK, Hubbard-Lucey VM, Sheehan KCF, Campbell KM, Lamb A, Ward JP, Sidney J, Blazquez AB, Rech AJ, Zaretsky JM, Comin-Anduix B, Ng AHC, Chour W, Yu TV, Rizvi1 H, Chen JM, Manning P, Steiner GM, Doan XC, The TESLA Consortium, Merghoub T, Guinney J, Kolom A, Selinsky C, Ribas A, Hellmann MD, Hacohen N, Sette A, Heath JR, Bhardwaj N, Ramsdell F, Schreiber RD, Schumacher TN, Kvistborg P, Defranoux N. Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction. Cell. 2020.

## Installation

Two methods exist to run `antigen.garnish`:

1. Docker
2. Linux

### Docker

```sh
docker pull leeprichman/antigen_garnish

cID=$(docker run -it -d leeprichman/antigen_garnish /bin/bash)
```

[Download](https://services.healthtech.dtu.dk/software.php) netMHC binaries (academic license): NetMHC 4.0, NetMHCpan 4.1b, NetMHCII 2.3, NetMHCIIpan 4.0.

Copy netMHC `tar.gz` files to the container and run the installation script:

```sh
docker cp netMHC-4.0a.Linux.tar.gz cID:/netMHC-4.0a.Linux.tar.gz
docker cp netMHCII-2.3.Linux.tar.gz cID:/netMHCII-2.3.Linux.tar.gz
docker cp netMHCIIpan-4.0.Linux.tar.gz cID:netMHCIIpan-4.0.Linux.tar.gz
docker cp netMHCIIpan-4.0.Linux.tar.gz cID:netMHCIIpan-4.0.Linux.tar.gz

docker exec $cID config_netMHC.sh
```

### Linux

#### Dependencies

- R &ge; 3.5.0
- [Bioconductor](https://www.bioconductor.org/install/) [Biostrings package](https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- [mhcflurry](https://github.com/openvax/mhcflurry)
- tcsh (required for netMHC, install via your package manager)

#### Installation

Install the dependencies listed above. Then, download and extract `antigen.garnish` data:

```sh
ANTIGEN_GARNISH_DIR="~/antigen.garnish"

cd ~
curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz
chmod 700 -R "$ANTIGEN_GARNISH_DIR"
```

Install antigen.garnish:

```r
# install.packages("devtools")
devtools::install_github("immune-health/antigen.garnish")
```

Next, [download](https://services.healthtech.dtu.dk/software.php) netMHC binaries (academic license): NetMHC 4.0, NetMHCpan 4.1b, NetMHCII 2.3, NetMHCIIpan 4.0.

Move the binaries into the `antigen.garnish` data directory, first setting the `NET_MHC_DIR` and `ANTIGEN_GARNISH_DIR` environment variables:

```sh
NET_MHC_DIR=/path/to/folder/containing/netMHC/downloads

cd "$NET_MHC_DIR"
mkdir -p "$ANTIGEN_GARNISH_DIR/netMHC"

tar xvzf netMHC-4.0a.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"
tar xvzf netMHCII-2.3.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"
tar xvzf netMHCIIpan-4.0.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"
tar xvzf netMHCIIpan-4.0.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"

chown "$USER" "$ANTIGEN_GARNISH_DIR/netMHC"
chmod 700 -R "$ANTIGEN_GARNISH_DIR/netMHC"
```

## Usage

See the [website](https://neoantigens.rech.io/reference/index.html) or [reference manual](https://get.rech.io/antigen.garnish.pdf).

### Docker

#### Interactive use

```sh

cID=$(docker run -it -d leeprichman/antigen_garnish /bin/bash)
docker exec -it $cID bash
R
```

```r
library("antigen.garnish")
```

#### VCF input

Paired tumor-normal VCFs annotated with SnpEff against any GRCh38 or GRCm38 releases are supported. For many variant callers, the tumor sample name is "TUMOR". In this case, the following input VCF file names will work:

```
TUMOR.vcf
TUMOR_se.vcf
TUMOR.ann.vcf
TUMOR.vcf.gz
```

MHC input is a JSON file from xHLA or a 2-column tab or comma-separated file ending in "mhc.txt" with the following format:

```
Example mouse .csv file:

		sample_id,MHC
		mysample.vcf,H-2-Kb H-2-Db H-2-IAb

Example human .csv file:

		sample_id,MHC
		mysample.vcf,HLA-A*02:01 HLA-B*07:02 or H-2-Kb H-2-Db
```

```sh
VCFO="TUMOR.vcf"
MHC="hla.json"

cID=$(docker run -it -d leeprichman/antigen_garnish /bin/bash)

docker cp $VCFO $cID:/$VCFO
docker cp $MHC $cID:/$MHC

# run antigen.garnish
docker exec $cID run_antigen.garnish.R
```

#### Peptide or transcript input

Start and configure the container as described above. One or more tab or comma-separated input files with the pattern "docker_input" in the file name with the following format are required:

```
Example transcript input .csv file:

		sample_id,ensembl_transcript_id,cDNA_change,MHC
		sample_1,ENST00000311936,c.718T>A,HLA-A*02:01 HLA-A*03:01

Example peptide input .csv file:

		sample_id,pep_mut,mutant_index,MHC
		sample_1,MTEYKLVVVDAGGVGKSALTIQLIQNHFV,25,HLA-A*02:01 HLA-A*03:01
```

```sh
INPUT="docker_input.csv"

cID=$(docker run -it -d leeprichman/antigen_garnish /bin/bash)

docker cp $INPUT $cID:/$DT

# run antigen.garnish
docker exec $cID run_antigen.garnish_direct.R
```

### Linux

#### Predict neoantigens from missense mutations, insertions, and deletions

```r
library(magrittr)
library(data.table)
library(antigen.garnish)

# load an example VCF
dir <- system.file(package = "antigen.garnish") %>%
       file.path(., "extdata/testdata")

file <- file.path(dir, "TUMOR.vcf")

# extract variants
dt <-  garnish_variants(file)

# add space separated MHC types
# see list_mhc() for nomenclature of supported alleles
# MHC may also be set to "all_human" or "all_mouse" to use all supported alleles

dt[, MHC := c("HLA-A*01:47 HLA-A*02:01 HLA-DRB1*14:67")]

# predict neoantigens
result <- dt %>% garnish_affinity(.)

result %>% str
```

#### Directly calculate foreignness score and dissimilarity for a list of sequences

```r
library(magrittr)
library(data.table)
library(antigen.garnish)

# generate our character vector of sequences
v <- c("SIINFEKL", "ILAKFLHWL", "GILGFVFTL")

# calculate foreignness score
v %>% foreignness_score(db = "human") %>% print

# calculate dissimilarity
v %>% dissimilarity_score(db = "human") %>% print
```

#### How are peptides generated?

```r
library(magrittr)
library(data.table)
library(antigen.garnish)

data.table::data.table(
   pep_base = "Y___*___THIS_IS_________*___A_PEPTIDE_TEST!______*__X",
   mutant_index = c(5, 25, 47, 50),
   pep_type = "test",
   var_uuid = c(
                "front_truncate",
                "middle",
                "back_truncate",
                "end")) %>%
   make_nmers %>% print
```

## Acknowledgments

We thank the follow individuals for contributions and helpful discussion:

- [David Balli](https://www.linkedin.com/in/davidballi1)
- [Adham Bear](https://www.med.upenn.edu/apps/faculty/index.php/g20001100/p1073)
- [Katelyn Byrne](https://www.parkerici.org/person/katelyn-byrne/)
- [Danny Wells](http://dannykwells.com/)

## License

Please see [LICENSE](https://github.com/immune-health/antigen.garnish/blob/master/LICENSE).
