## Detailed install instructions.

Recommended instructions for installation of `antigen.garnish`.  Install onto an AWS EC2 instance running an Ubuntu Server 16.04 LTS per-configured with [Bioconductor](https://www.bioconductor.org/help/bioconductor-cloud-ami/#overview) (ami-aab1e9d0) and the top 80 bioconductor [tools](http://www.bioconductor.org/packages/stats/).  This AMI runs R version 4.3.2 and Bioconductor version 3.6.

#### 1. Update ubuntu

```
sudo apt-get update
```

##### 6. Install [antigen.garnish](https://github.com/andrewrech/antigen.garnish) and [dt.inflix](https://github.com/andrewrech/dt.inflix) in R and run tests.

```
install.packages(c("devtools", "testhat"))

library("devtools")
library("testthat")

devtools::install_github(c("andrewrech/dt.inflix", "andrewrech/antigen.garnish", "tidyverse/magrittr"), force = TRUE)
```

##### 2. Download `antigen.garnish` dependency files
```
cd ~
curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz
chmod 777 -R ./antigen.garnish
chown `whoami` ./antigen.garnish
sudo mv ./antigen.garnish/ncbi-blast-2.7.1+/bin/* /usr/local/bin
```

#### 3. Install [mhcflurry](https://github.com/openvax/mhcflurry) and download mhcflurry prediction models.

```
sudo pip --disable-pip-version-check install mhcflurry
mhcflurry-downloads fetch
```

#### 4. Test `antigen.garnish` Installation in R.

```
testthat::test_package("antigen.garnish")
```
