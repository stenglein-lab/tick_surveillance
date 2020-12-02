#!/bin/bash 

# This script runs the commands necessary to setup a linux server to run
# the CDC Tick Surveillance Pipeline
#
# Mark Stenglein
# 12/1/2020

cd $HOME
curl -OL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x $HOME/Miniconda3-latest-Linux-x86_64.sh 

echo "The miniconda setup program will now run to install miniconda."
echo " "
echo "miniconda is a minimal conda installation.  Conda is a system for installing and"
echo "managing software dependencies."  
echo " "
echo "For more information, see: https://docs.conda.io/en/latest/miniconda.html#linux-installers"
echo " "
echo "Enter to continue."
read x

# run the miniconda setup script
$HOME/Miniconda3-latest-Linux-x86_64.sh

# conda init
conda init

# install nextflow in base conda environment
conda install -c bioconda nextflow=20.10.*

# create a conda environment
conda env create --prefix $HOME/tick_conda_environment -f tick_conda_environment.yaml

# use shorter conda environment names 
# see https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#specifying-a-location-for-an-environment
conda config --set env_prompt '({name})'

# install some R packages into the new conda environment
conda activate $HOME/tick_conda_environment

# install some system packages necessary for the subsequent R packages:
# sudo apt-get install libcurl4-openssl-dev
# sudo apt-get install libxml2-dev

# install some necessary R packages
# Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'
# Rscript -e 'install.packages("dada2", repos="https://cloud.r-project.org")'
# Rscript -e 'install.packages("openxlsx", repos="https://cloud.r-project.org")'

