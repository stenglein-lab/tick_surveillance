#!/usr/bin/env Rscript 
#
# This script runs dada2 on a set of trimmed fastq files 
# 
# See: 
# https://benjjneb.github.io/dada2/tutorial.html
#

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
  trimmed_path=args[1]
  outdir="./"
} else {
  # if running via RStudio (for development or troubleshooting)
  r_bindir = "."
  trimmed_path = NA
  outdir="../results/"
}

# get lists of fastq files in trimmed directory 
fnFs <- sort(list.files(trimmed_path, pattern="_R1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(trimmed_path, pattern="_R2_trimmed.fastq.gz", full.names = TRUE))

# extract sample names from fastq file names
# so that sample names will match execpted values provided to sample sheet 
# when demultiplexing

# remove _R[12]_trimmed.fastq.gz from file names to make better sample names
sample.names <- gsub("_R[12]_trimmed.fastq.gz", "", basename(fnFs))


# This will place dada-filtered files in dada_filtered/ subdirectory of the pwd
filtFs <- file.path("dada_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("dada_filtered", paste0(sample.names, "_R_filt.fastq.gz"))


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


# only keep datasets with >0 reads to avoid dada2 throwing an error because of empty datasets
# empty datasets will not contribute to error learning and will have no assigned ASVs, as expected.
# see: https://github.com/benjjneb/dada2/issues/469
keep <- out[,"reads.out"] > 10 
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]


# also only keep file pairs for which the filtered files exist
# see: https://github.com/benjjneb/dada2/issues/375
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

# Could save this as output...
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# save the error profile for bases with a Q-score of 35
# q35_errors <- errF$err_out[,36]
# q35_errors
# write.table(q35_errors, "q35_errors.txt", sep="\t")

# Run dada function to identify ASVs
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# create a table of dada2 results
seqtab <- makeSequenceTable(mergers)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


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

# write the output to a table
write.table(t, paste0(outdir, "dada_seqtab.txt"), sep="\t", col.names=T, quote=F)

# write out version info into versions.yml
writeLines(c("DADA2:", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")),paste0("    ShortRead: ", packageVersion("ShortRead")) ), "versions.yml")


# Track reads and create output file (https://benjjneb.github.io/dada2/tutorial.html)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# write read tracking info to csv file
write.csv(track, paste0(outdir, "dada_read_clean_all.csv"), col.names=T)
