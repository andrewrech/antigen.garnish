## Detailed install instructions.

Recommended instructions for installation of `antigen.garnish`.  Install onto an AWS EC2 instance running an Ubuntu Server 16.04 LTS per-configured with [Bioconductor](https://www.bioconductor.org/help/bioconductor-cloud-ami/#overview) (ami-aab1e9d0) and the top 80 bioconductor [tools](http://www.bioconductor.org/packages/stats/).  This AMI runs R version 4.3.2 and Bioconductor version 3.6.

#### Requirements for installation

* Linux
* R ≥ 3.4
* python2.7 and pip

#### 1. Update ubuntu
```
sudo apt-get update
```

##### 2. Install [antigen.garnish](https://github.com/andrewrech/antigen.garnish) and [dt.inflix](https://github.com/andrewrech/dt.inflix) in R.
```
Rscript --vanilla -e \
'install.packages("devtools", repos = "http://cran.us.r-project.org"); devtools::install_github("hadley/devtools"); install.packages("testthat", repos = "http://cran.us.r-project.org")'

Rscript --vanilla -e \
'devtools::install_github(c("tidyverse/magrittr", "andrewrech/dt.inflix", "andrewrech/antigen.garnish"))'
```

##### 3. Download `antigen.garnish` dependency files.  

* Will download peptide and cDNA databases, known immunogenic IEDB sequences, [NCBI blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [mhcnuggers](https://www.biorxiv.org/content/biorxiv/early/2017/06/23/154757.full.pdf), and [netMHC](http://www.cbs.dtu.dk/services/software.php) tools.  Please note that netMHC, netMHCpan, netMHCII, and netMHCIIpan require academic-use only licenses.

```
cd ~
curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz
chmod 777 -R ./antigen.garnish
chown `whoami` ./antigen.garnish
sudo mv ./antigen.garnish/ncbi-blast-2.7.1+/bin/* /usr/local/bin
```

#### 4. Install [mhcflurry](https://github.com/openvax/mhcflurry) and download mhcflurry prediction models.
```
sudo pip --disable-pip-version-check install scipy mhcflurry h5py biopython

mhcflurry-downloads fetch
```

#### 5. Test `antigen.garnish` Installation in R.
```r
testthat::test_package("antigen.garnish")
```
