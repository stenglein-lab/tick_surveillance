#
# This script assigns observed sequences to a set of reference sequences based
# on blast alignments
#

library(tidyverse)

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
  tidy_table_path=args[2]
  blast_output_path=args[3]
} else {
  # if running via RStudio
  r_bindir = "."
  tidy_table_path="../results/tidy_sequence_table.tsv"
  blast_output_path="../results/observed_sequences.fasta.bn_refseq"
}

# read in the tidy_sequecne_table
tidy_sequence_table <- read.delim(tidy_table_path, sep="\t", header=T)

# make a dataframe containing the sequences and their sequence number (which matches the blast query)
# and add a column listing their length
sequences <- tidy_sequence_table %>% 
  group_by(sequence, sequence_number) %>% 
  summarize(.groups = "drop") %>%
  mutate(sequence_length = str_length(sequence))
  

blast_df <- read.delim(blast_output_path, sep="\t", header=F)
colnames(blast_df) <- c("query", "subject", "percent_identity", "alignment_length", "mismatches", "gapopens", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# consolidated_blast_df_all <- blast_df %>% group_by(query) %>% arrange(-bitscore, .by_group = T) 

# take the highest scoring blast hit for each query
consolidated_blast_df <- blast_df %>% group_by(query) %>% arrange(-bitscore, .by_group = T) %>% filter(row_number() == 1)

# merge the sequence info into the blast output
consolidated_blast_df <- left_join(consolidated_blast_df, sequences, by=c("query" = "sequence_number"))

# what fraction of the query was present in the alignment?
consolidated_blast_df <- consolidated_blast_df %>% mutate(fraction_query_aligned = 100 * alignment_length / sequence_length)

ggplot(consolidated_blast_df) + geom_point(aes(x=fraction_query_aligned, y=percent_identity))

# --------------------------------------------------------------------------
# logic for saying that one of the observed sequences is close enough 
# to one of the reference sequences to assign it to that reference sequence
# --------------------------------------------------------------------------

min_percent_identity <- 96
min_fraction_aligned <- 90

assigned_blast_df <- consolidated_blast_df %>% filter(percent_identity > min_percent_identity & fraction_query_aligned > min_fraction_aligned)

assigned_queries <- assigned_blast_df %>% select(query)

unassigned_sequences <- is_in(sequences$sequence_number, assigned_queries$query)

sequences[!(sequences$sequence_number %in% assigned_queries$query),]
