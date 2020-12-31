![](https://s3.amazonaws.com/get.rech.io/antigen.garnish_build_status.svg) [![rech.io](https://img.shields.io/badge/endpoint.svg?url=https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.json)](https://s3.amazonaws.com/get.rech.io/antigen.garnish_coverage.html) ![](https://img.shields.io/badge/docker-andrewrech/antigen.garnish:2.1.1-brightgreen)

# antigen.garnish 2

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
docker pull andrewrech/antigen.garnish

cID=$(docker run -it -d andrewrech/antigen.garnish /bin/bash)
```

[Download](https://services.healthtech.dtu.dk/software.php) netMHC binaries (academic license): NetMHC 4.0, NetMHCpan 4.1b, NetMHCII 2.3, NetMHCIIpan 4.0.

Copy netMHC `tar.gz` files to the container and run the installation script:

```sh
docker cp netMHC-4.0a.Linux.tar.gz $cID:/netMHC-4.0a.Linux.tar.gz
docker cp netMHCII-2.3.Linux.tar.gz $cID:/netMHCII-2.3.Linux.tar.gz
docker cp netMHCpan-4.1b.Linux.tar.gz $cID:netMHCpan-4.1b.Linux.tar.gz
docker cp netMHCIIpan-4.0.Linux.tar.gz $cID:netMHCIIpan-4.0.Linux.tar.gz

docker exec $cID config_netMHC.sh
```

### Linux

#### Dependencies

- R &ge; 3.5.0
- [Bioconductor](https://www.bioconductor.org/install/) [Biostrings package](https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- [mhcflurry](https://github.com/openvax/mhcflurry)
- netMHC tools (see below)
- tcsh (required for netMHC, install via your package manager)

#### Installation

Install the dependencies listed above. Then, download and extract `antigen.garnish` data:

```sh
ANTIGEN_GARNISH_DIR="~/antigen.garnish"

cd ~
curl -fsSL "http://get.rech.io/antigen.garnish-2.1.0.tar.gz" | tar -xvz
chmod 700 -R "$ANTIGEN_GARNISH_DIR"
```

Install antigen.garnish:

```r
# install.packages("remotes")
remotes::install_github("andrewrech/antigen.garnish")
```

Next, [download](https://services.healthtech.dtu.dk/software.php) netMHC binaries (academic license): NetMHC 4.0, NetMHCpan 4.1b, NetMHCII 2.3, NetMHCIIpan 4.0.

Move the binaries into the `antigen.garnish` data directory, first setting the `NET_MHC_DIR` and `ANTIGEN_GARNISH_DIR` environment variables:

```sh
NET_MHC_DIR=/path/to/folder/containing/netMHC/downloads

cd "$NET_MHC_DIR"
mkdir -p "$ANTIGEN_GARNISH_DIR/netMHC"

tar xvzf netMHC-4.0a.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"
tar xvzf netMHCII-2.3.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"
tar xvzf netMHCpan-4.1b.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"
tar xvzf netMHCIIpan-4.0.Linux.tar.gz -C "$ANTIGEN_GARNISH_DIR/netMHC"

chown "$USER" "$ANTIGEN_GARNISH_DIR/netMHC"
chmod 700 -R "$ANTIGEN_GARNISH_DIR/netMHC"
```

## Usage

See the [reference manual](https://get.rech.io/antigen.garnish.pdf).

### Docker

#### Interactive use

Copy any VCF files and/or metadata including HLA alleles onto the running container using the `docker cp` [command](https://docs.docker.com/engine/reference/commandline/cp/). The container ID is still saved as `$cID` from the installation above. You will also need to use this command and container ID to copy saved output files from the docker container after you complete your analysis.

Copy any needed files onto the running container, for example:

```sh

docker cp myfile.txt $cID:/myfilecopy.txt

```

Now launch the interactive virtual machine with the container you started:

```sh

docker exec -it $cID bash
R

```

```r

library(antigen.garnish)

```

Follow the instructions in the next section titled **Linux** to complete your interactive R analysis. When you complete your analysis, copy any desired output files off the container to your local machine with the `docker cp` [command](https://docs.docker.com/engine/reference/commandline/cp/). Shut down and clean up your container like this:

```sh

docker cp $cID:/myoutput.txt ~/myagdockeroutput.txt

docker stop $cID

docker rm $cID

```

### Linux

Parallel cores used can be set via environment variable AG_THREADS (default: all available).

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

#### Predict neoantigens from Microsoft Excel or other table input

Transcript ID level input table format:
```
# sample_id ensembl_transcript_id cDNA_change MHC
# sample_1  ENST00000311936       c.718T>A    HLA-A*02:01 HLA-A*03:01
# sample_1  ENST00000311936       c.718T>A    H-2-Kb H-2-Kb
```

Protein level input (with optional WT paired input) table format:
```
# sample_id pep_mut           pep_wt            mutant_index MHC
# sample_1  MTEYKLVVVDADGVGK  MTEYKLVVVDAGGVGK  12           HLA-A*02:01
# sample_1  MTEYKLVVVDDDGVGK  MTEYKLVVVDAGGVGK  12 13        HLA-A*02:01
# sample_1  MTEYKLVVVDAGGAAA  MTEYKLVVVDAGGVGK  14 15 16     HLA-A*02:01
# sample_1  SIINFEKLMILKATFI  MTEYKLVVVDAGGVGK  all          HLA-A*02:01
```

```r
library(magrittr)
library(data.table)
library(antigen.garnish)
library(rio) # package to import Excel and other tables

# load an example table
dir <- system.file(package = "antigen.garnish") %>%
       file.path(., "extdata/testdata")

file <- file.path(dir, "antigen.garnish_example_peptide_with_WT_input.txt")

# read in excel or other format file with rio::import and convert to data table
# or substitute the path to your file here
mytable <- rio::import(file) %>% data.table::as.data.table()

# only use first two rows of table for example
mytable <- mytable[1:2]

# predict neoantigens from data table object
result <- garnish_affinity(mytable)

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
