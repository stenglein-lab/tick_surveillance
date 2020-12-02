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
  trimmed_path = "../trimmed_fastq"
}

# get lists of fastq files in trimmed directory 
fnFs <- sort(list.files(trimmed_path, pattern="_R1_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(trimmed_path, pattern="_R2_trimmed.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

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

# Could save this as output...
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

# Ru 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# ?mergePairs

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
head(track)

# TODO: output this table for tracking purposes

# convert to a data frame
t <- as.data.frame(seqtab.nochim)

# make a new column based on row names
t$dataset <- rownames(t)

# this is now a tidy formatted table with:
# unique, merged sequences
# their abundances in each dataset
tidy_sequence_table <- pivot_longer(t, -dataset, names_to = "sequence", values_to = "abundance" )

# pivot wider
wide_sequence_table <- pivot_wider(tidy_sequence_table, names_from = sequence, values_from = abundance)
write.table(wide_sequence_table, "wide_seqtab.txt", sep="\t", quote=F, row.names = F)

sequence_df <- data.frame(seq_id=integer(),
                 sequence=character(), 
                 stringsAsFactors=FALSE) 
# all the sequences in a vector
sequences <- distinct(tidy_sequence_table, sequence)
sequences <- sequences %>% mutate(sequence_number = row_number())
# write.table(sequences, "sequences.txt", sep="\t", quote=F, row.names = F)

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

writeFasta(sequences, "sequences.fasta")

# 9/30/2020
# TODO: 
# - assign sequences to one of the reference sequences
 # - what tool to use for this?   BLAST?   other?
 # - where to do this?  in R?  in a bash script or something?
 # - how to handle results?  
# - do something with non-matching sequences...




