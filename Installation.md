## Detailed install instructions.


Recommended instructions for installation of `antigen.garnish`. Install onto an AWS EC2 instance running an Ubuntu Server 16.04 LTS per-configured  with [Bioconductor](https://www.bioconductor.org/help/bioconductor-cloud-ami/#overview) and the top 80 [tools](http://www.bioconductor.org/packages/stats/).  This AMI runs R version 4.3.2 and Bioconductor version 3.6


#### 1. Update and install build essentials.


```
sudo apt-get update
sudo apt-get install build-essential, zlib1g-dev  
```

#### 2. Download and install [conda](https://repo.continuum.io/).  Will require license agreement and setting environment path for anaconda directory.


```
wget https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh
bash Anaconda2-5.1.0-Linux-x86_64.sh
# after successful installation
source ~/.bashrc
```

#### 3. Install dependencies from conda.

```
conda update -n base conda
conda install -c r r-testthat 
conda install -c r r-devtools
conda install -c anaconda libcurl 
conda install -c conda-forge scipy 
conda install -c conda-forge h5py

```

#### 4. Install [mhcflurry](https://github.com/openvax/mhcflurry) and download mhcflurry prediction models. 


```
pip --disable-pip-version-check install mhcflurry 
mhcflurry-downloads fetch
```

##### 5. Download `antigen.garnish` dependency files 

```
cd ~
curl -fsSL "http://get.rech.io/antigen.garnish.tar.gz" | tar -xvz
chmod 777 -R ./antigen.garnish
chown `whoami` ./antigen.garnish
sudo mv ./antigen.garnish/ncbi-blast-2.7.1+/bin/* /usr/local/bin
```

##### 6. Install [antigen.garnish](https://github.com/andrewrech/antigen.garnish) and [dt.inflix](https://github.com/andrewrech/dt.inflix) in R and run tests.

```
library("devtools")
library("testthat")

devtools::install_github(c("andrewrech/dt.inflix", "andrewrech/antigen.garnish"))

testthat::test_package("antigen.garnish")
```

