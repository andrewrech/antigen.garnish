# antigen.garnish on AWS

Recommended instructions for minimal installation of `antigen.garnish` onto a fresh AWS EC2 instance running Ubuntu Server 18.04 LTS (ami-0ac019f4fcb7cb7e6). This image has many of the necessary dependencies installed, including python (for MHCflurry).

## Bootstrap an AWS instance

From the instance command line:

#### 1. add the apt CRAN repository for R

```sh
cd "$HOME"

# add the R Ubuntu bionic apt repository
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/'

# update ubuntu with the added repo
sudo apt-get update -y
```

#### 2. install R, pip package manager, and miscellaneous dependencies/libraries

```sh
sudo apt-get install -y r-base \
 python-pip \
 libcurl4-gnutls-dev \
 libssl-dev subversion \
 libxml2-dev \
 htop \
 libbz2-dev \
 liblzma-dev
```
#### 3. run antigen.garnish install script
This shell script, also available in the [github repository](https://github.com/immune-health/antigen.garnish/blob/master/inst/extdata/install.sh), wraps and installs the MHC prediction tools, additional R packages, including from [Bioconductor](https://www.bioconductor.org/), and downloads the antigen.garnish databases of human and murine transcript metadata, known immunogenic sequences, and self-sequences.

```sh
curl -fsSL http://get.rech.io/antigen.garnish.sh | sudo sh
```

#### 4. test `antigen.garnish`
You can now test antigen.garnish, run the below from the command line. Once this returns no failures, you may use antigen.garnish from R interactively or from the shell using Rscripts.

```r
sudo Rscript -e 'testthat::test_package("antigen.garnish")'
```

## Still having trouble?
If this isn't working with your server or you do not have the sudo permissions necessary to install, consider running antigen.garnish with [Docker](https://www.docker.com/get-started). Multiple docker containers can be orchestrated to analyze a large number of samples in parallel in a 100% sealed and consistent environment. Detailed instructions for running antigen.garnish with Docker can be found [here](https://github.com/immune-health/antigen.garnish/wiki/Docker).
