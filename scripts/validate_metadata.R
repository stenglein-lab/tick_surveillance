#
# This script validates the metadata file  
#
# Mark Stenglein 2/24/2022
#

library(tidyverse)
library(readxl)

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#

# if running from Rscript (called from the pipeline)
if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  r_bindir = args[1]
  sample_metadata_file = args[2]
  sample_ids_file = args[3]
  output_dir = "./"
} else {
  # if running via RStudio
  r_bindir  =  "."
  sample_metadata_file = "../input/AK_metadata.xlsx"
  sample_ids_file = "../results/sample_ids.txt"
  output_dir = "../results/"
}


# aspects of metadata to validate:
# 
# 1. that the appropriate columns are included
# 2. that sample IDs of fastq files match one sample ID in the metadata file
# 3. warn if there are sample IDs in the metadata that don't have corresponding fastq files 


