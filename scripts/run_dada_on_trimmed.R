#
# This script runs dada2 on a set of trimmed fastq files 
# 
# See: 
# https://benjjneb.github.io/dada2/tutorial.html
#

library(tidyverse)
library(dada2)

#
# This code block sets up input arguments to either come from the command line
# (if not running in interactive mode) or to use expected default values 
# (if running in interactive mode, i.e. RStudio).  
#
# This is useful in case you want to modify this script
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  r_bindir=args[1]
  trimmed_path=args[2]
} else {
  # if running via RStudio
  r_bindir = "."
  trimmed_path = "../results/trimmed_fastq"
}

# get lists of fastq files in trimmed directory 
fnFs <- sort(list.files(trimmed_path, pattern="_R1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(trimmed_path, pattern="_R2_trimmed.fastq.gz", full.names = TRUE))

# extract sample names from fastq file names
# so that sample names will match execpted values provided to sample sheet 
# when demultiplexing

# remove _R[12]_trimmed.fastq.gz from file names to make better sample names
sample.names <- str_replace(basename(fnFs), "_R[12]_trimmed.fastq.gz", "")
# remove _S1_L001  etc. from file names to make better sample names
sample.names <- str_replace(sample.names, "_S\\S+_L[01]{3}", "")


# This will place dada-filtered files in dada_filtered/ subdirectory
filtFs <- file.path(trimmed_path, "dada_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(trimmed_path, "dada_filtered", paste0(sample.names, "_R_filt.fastq.gz"))


names(filtFs) <- sample.names
names(filtRs) <- sample.names

#
# Run dada2 filtering and trimming

# filterAndTrim options:
#
# truncQ 	
# (Optional). Default 2. Truncate reads at the first instance of a quality score less 
# than or equal to truncQ.

# truncLen 	
# (Optional). Default 0 (no truncation). Truncate reads after truncLen bases. 
# Reads shorter than this are discarded.

# maxEE
# (Optional). Default Inf (no EE filtering). After truncation, reads with higher than 
# maxEE "expected errors" will be discarded. Expected errors are calculated from the 
# nominal definition of the quality score: EE = sum(10^(-Q/10))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


# only keep file pairs for which the filtered files exist
# see: https://github.com/benjjneb/dada2/issues/375
# exists is a logical vector based on whether the R1 and R2 files exist
exists <- file.exists(filtFs) & file.exists(filtRs)
# use exists logical vector to subset filtered file names
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

# Estimate error rates: actual error rate as a function of quality score for each type   
# of base substitution.
#
# From the dada2 documentation (https://benjjneb.github.io/dada2/tutorial.html): 
#
# "The DADA2 algorithm makes use of a parametric error model (err) and every amplicon 
# dataset has a different set of error rates. The learnErrors method learns this error 
# model from the data, by alternating estimation of the error rates and inference of 
# sample composition until they converge on a jointly consistent solution."
#
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

errF$err_out

# Could save this as output...
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# save the error profile for bases with a Q-score of 35
# q35_errors <- errF$err_out[,36]
# q35_errors
# write.table(q35_errors, "q35_errors.txt", sep="\t")



# Ru 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
# dim(seqtab)

# table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab)

# write.table(seqtab, "65_seqtab.txt", sep="\t")

getN <- function(x) sum(getUniques(x))

# this will work for one sample
# getN(dadaFs)
# getN(dadaRs)
# getN(mergers)
# getN(seqtab.nochim)

# this works when >1 sample
# TODO: leave for most cases...
# sum(getUniques(dadaFs))
sapply(dadaFs, getN)
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
# head(track)

# COULD_DO: output the track table for tracking purposes : how many reads were collapsed/filtered at each step by dada

#
# Handle empty datasets.  These should be reported as datasets with abundances of 0 for all sequences.
#
# These could correspond to datasets that were truly empty after sequencing or had no reads after filtering.
#
# Add rows with all 0s for datasets that had no reads after filtering 
# this corresponds to the FALSE values in the exists vector
empty_sample_names <- sample.names[!exists]
# create a matrix with all 0s
new_rows_with_zeros <- matrix(0L, nrow = length(empty_sample_names), ncol = ncol(seqtab.nochim)) 
# name them according to their original names
row.names(new_rows_with_zeros) <- empty_sample_names
# concatenate this new all-zeros table onto the end of the data seqtab
seqtab.nochim <- rbind(seqtab.nochim, new_rows_with_zeros)

# convert to a data frame
t <- as.data.frame(seqtab.nochim)

# make a new column based on row names
t$dataset <- rownames(t)

# this is now a tidy formatted table with:
# unique, merged sequences
# and their abundances in each dataset
tidy_sequence_table <- pivot_longer(t, -dataset, names_to = "sequence", values_to = "abundance" )

# create a unique sequence id (#) for each sequence: for keeping track of in further processing steps
sequences <- tidy_sequence_table %>% group_by(sequence) %>% summarize () %>% mutate(sequence_number = row_number())

# join these sequence #s back into the tidy table
tidy_sequence_table <- left_join(tidy_sequence_table, sequences, by="sequence")

# write out the tidy-formatted table
write.table(tidy_sequence_table, "sequence_abundance_table.tsv", sep="\t", col.names=T, row.names=F)

# this writes fasta from the sequences data frame 
# see: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"sequence_number"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# write out the sequences in fasta format
writeFasta(sequences, "observed_sequences.fasta")





