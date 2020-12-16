#!/bin/bash 

# This script runs the commands necessary to setup a linux server to run
# the CDC Tick Surveillance Pipeline
#
# Mark Stenglein
#
# 12/1/2020

# current working directory
present_working_directory=`pwd`

# download Miniconda setup script

# change to home directory: will install the environment there
cd $HOME
curl -OL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x $HOME/Miniconda3-latest-Linux-x86_64.sh 

echo "The miniconda setup program will now run to install miniconda."
echo " "
echo "miniconda is a minimal conda installation.  Conda is a system for installing and"
echo "managing software dependencies."  
echo " "
echo "For more information, see: "
echo " https://docs.conda.io/en/latest/miniconda.html#linux-installers"
echo " "
echo "Enter to continue."
read x

# run the miniconda setup script
echo "running the miniconda installation script.  This may take several minutes to complete."
$HOME/Miniconda3-latest-Linux-x86_64.sh

# conda init
echo "initializing conda"
conda init

# install nextflow in base conda environment
echo "installing the nextflow pacakge in the base conda environment."
conda install -c bioconda nextflow=20.10.*

# create a conda environment
echo "creating a new conda environment with the software dependenceies"
echo "necessary to run the tick surveillance pipeline."
echo " "
echo "this will run for several minutes."
echo " "
echo "enter to continue"
read x

conda env create --prefix $HOME/tick_conda_environment -f ${present_working_directory}/tick_conda_environment.yaml

# use shorter conda environment names 
# see https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#specifying-a-location-for-an-environment
conda config --set env_prompt '({name})'

echo "conda environment setup complete."

# activate this conda enviornment so can use update_blastdb.pl
conda activate $HOME/tick_conda_environment

# download the NCBI nt database
echo "Downloading the NCBI nt database.  This will likely take several hours."
echo "enter to continue"
read x

update_blastdb.pl --decompress --num_threads 12 nt 
mkdir -p $HOME/nt_database
cp -Rp nt.* $HOME/nt_database





