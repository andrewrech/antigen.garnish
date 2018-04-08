## Detailed install instructions.


Instructions for installation of `antigen.garnish` on an AWS EC2 instance running an Ubuntu Server 16.04 LTS (ami-43a15f3e).  


1. Update and install build essentials.


`
sudo apt-get update
sudo apt-get install build-essential, zlib1g-dev  
`

2. Download and install conda.  Will require license agreement and setting environment path for anaconda directory.


`
wget https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh
bash Anaconda2-5.1.0-Linux-x86_64.sh

# after successful installation
source ~/.bashrc
`

3. Install dependencies from conda.


`
conda update -n base conda
conda install -c r r-essentials
conda install -c r r-testthat 
conda install -c r r-devtools
conda install -c r r-git2r
conda install -c conda-forge scipy 
conda install -c conda-forge h5py
`

4. Install [mhcflurry](https://github.com/openvax/mhcflurry) and download mhcflurry prediction models. 


`
pip --disable-pip-version-check install mhcflurry 
mhcflurry-downloads fetch
`
