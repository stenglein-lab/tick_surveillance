#!/bin/bash

# download and unpack the NCBI nt database
# Mark Stenglein 1.9.2021

# Doing it this way because update_blastdb.pl [from NCBI] continually dropped connection

# create a list of nt files
# assume the dbs will go up to 99 sub-files - may have to up this in the future
# if goes beyond 99 will have to make sure numbering is ok, w/ regards to 0 padding...

# have to do this cd command in case this script is run as a cron job
# in which case it would normally excute w/ home dir as pwd
DATABASE_DIR=$HOME/NCBI_databases/nt/
mkdir -p $DATABASE_DIR
cd $DATABASE_DIR

echo "***************************************************************"
echo "** downloading taxonomy info so blast will be taxonomy aware **
echo "***************************************************************"


# first, get and unpack the taxdb file blast cli applications need to be 'taxonomically aware'
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar xvf taxdb.tar.gz
rm -f taxdb.tar.gz

# setup BLASTDB environmental variable so blast knows where to find nt database and taxonomy info
# delete any such matching lines
sed -i '/BLASTDB/d' ~/.bashrc
echo "export BLASTDB=$DATABASE_DIR" >> ~/.bashrc

echo "***************************************************"
echo "** downloading and updating nt files and indexes **"
echo "***************************************************"
date

num_subfiles=99

subfile_nums=(`seq -s " " -w 0 $num_subfiles`)

nt_files=()

subfile=1
for subfile in ${subfile_nums[@]}
do
   nt_file="nt.${subfile}.tar.gz"
   nt_files+=($nt_file)
done

echo "starting download and unpacking of files"
date

echo "fetching nt files from NCBI FTP server"
for f in ${nt_files[@]}
do
   wget -N ftp://ftp.ncbi.nlm.nih.gov/blast/db/$f
done

echo "unpacking nt files"
for f in `ls nt.??.tar.gz`
do
   tar xfzv $f
done

echo "finished downloading and unpacking files"
date

date >> nr_nt_updated_dates.txt

# remove tar.gz files
echo "deleting .tar.gz files"
rm -f nt*.tar.gz


