FROM ubuntu:20.04 as dependencies

# to run tests
# export DOCKER_BUILDKIT=1 BUILDKIT_PROGRESS=plain; docker build --build-arg CACHEBUST="$(date +%s)" --target test -t andrewrech/antigen.garnish -f Dockerfile .

# build documentation
# export DOCKER_BUILDKIT=1 BUILDKIT_PROGRESS=plain; docker build --build-arg CACHEBUST="$(date +%s)" --target docs -t andrewrech/antigen.garnish -f Dockerfile .

ENV DEBIAN_FRONTEND=noninteractive
ENV TERM linux
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        default-jre \
        git \
        gnupg \
        libarchive-tools \
        libbz2-dev \
        libcurl4-gnutls-dev \
        libidn11 \
        liblzma-dev \
        liblzo2-2 \
        libssl-dev \
        libxml2-dev \
        locales \
        parallel \
        python-setuptools \
        python3-pip \
        r-base \
        r-base-dev \
        r-recommended \
        software-properties-common \
        subversion \
        tcsh \
    && rm -rf /var/lib/apt/lists/*

## https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_BASE_VERSION 3.5.2
ENV DEBIAN_FRONTEND=noninteractive
ENV TERM linux

RUN apt-key adv --keyserver "hkp://keyserver.ubuntu.com:80" \
                --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository \
       'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu bionic-cran40/'

# install python dependencies
RUN pip3 --disable-pip-version-check \
      install --no-cache-dir \
      scipy==1.4.1 \
      h5py \
      mhcflurry \
      biopython \
    && mhcflurry-downloads fetch

## mhcflurry does not work as root user by default
## copy needed files to prevent errors
RUN mkdir /.local/ &&  \
      mkdir /.local/share && \
      mv /root/.local/share/mhcflurry /.local/share/

RUN Rscript --vanilla -e \
    'install.packages(c("BiocManager", "testthat", "rcmdcheck", "data.table", "mclust", "Rdpack", "roxygen2", "tidyr", "uuid", "vcfR", "zoo"), repos = "http://cran.us.r-project.org"); BiocManager::install("Biostrings")'

FROM dependencies as data
ARG ANTIGEN_GARNISH_DATA_LINK=https://get.rech.io/antigen.garnish-2.1.0.tar.gz
# ARG ANTIGEN_GARNISH_DATA_LINK=https://get.rech.io/antigen.garnish-dev-no-distrib.tar.gz
ARG CACHEBUST_DATA
WORKDIR /root
RUN mkdir -p ./antigen.garnish \
       && curl -fsSL "$ANTIGEN_GARNISH_DATA_LINK" > antigen.garnish.tar.gz \
       && tar xvf antigen.garnish.tar.gz -C ./antigen.garnish

FROM dependencies as install
ARG CACHEBUST_DATA
ARG CACHEBUST
COPY --from=data /root/antigen.garnish/* /root/antigen.garnish/
WORKDIR /root/src
COPY . ./
RUN cp ./inst/extdata/src/config_netMHC.sh \
       /usr/local/bin \
       && R CMD INSTALL . \
       && echo 'export AG_DATA_DIR="/root/antigen.garnish"' >> /root/.bashrc \
    && echo 'export PATH=/root/antigen.garnish/netMHC/netMHC-4.0:/root/antigen.garnish/netMHC/netMHCII-2.3:/root/antigen.garnish/netMHC/netMHCIIpan-4.0/:/root/antigen.garnish/netMHC/netMHCpan-4.1:/root/antigen.garnish/ncbi-blast-2.10.1+-src/c++/ReleaseMT/bin:"$PATH"' >> /root/.bashrc

# stage for testing, skipped in production image
FROM install as test
WORKDIR /root/src
RUN Rscript -e 'pkg <- c("BiocManager", "testthat", "rcmdcheck", "data.table", "mclust", "Rdpack", "roxygen2", "tidyr", "uuid", "vcfR", "zoo"); installed <- pkg %in% installed.packages()[,1]; if (!all(installed)){print("Failed to install:"); print(pkg[!pkg %in% installed.packages()[,1]]); quit(save = "no", status = 1, runLast = FALSE)}'
RUN Rscript --vanilla -e 'source("/root/src/tests/testthat/setup.R"); testthat::test_dir("/root/src/tests/testthat", stop_on_failure = TRUE)'

# stage for documentation, skipped in production image
FROM install as docs
WORKDIR /root/src
RUN apt-get install -y --no-install-recommends \
      texinfo \
      && Rscript --vanilla -e 'install.packages("tinytex"), repos = "http://cran.us.r-project.org"); tinytex::install_tinytex()' \
      && export PATH=/root/bin:"$PATH" \
      && R CMD Rd2pdf --no-preview -o antigen.garnish.pdf .

FROM install as release
WORKDIR /root
CMD ["bash"]
