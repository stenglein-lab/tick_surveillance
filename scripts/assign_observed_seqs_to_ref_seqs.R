#
# This script assigns observed sequences to a set of reference sequences based
# on blast alignments
#
# Mark Stenglein 12/11/2020
#

library(tidyverse)
library(openxlsx)
library(readxl)
library(DT)

#

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  r_bindir=args[1]
  tidy_table_path=args[2]
  blast_output_path=args[3]
  sample_metadata_file=args[4]
  dataset_ids=args[5]
} else {
  # if running via RStudio
  r_bindir = "."
  tidy_table_path="../results/sequence_abundance_table.tsv"
  blast_output_path="../results/observed_sequences.fasta.bn_refseq"
  sample_metadata_file="../input/sample_metadata.xlsx"
}


# write out a file of unassigned sequence# this writes fasta from the sequences data frame 
# see: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
# TODO: import (source) 
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

# read in the sequence_abundance_table.tsv file, which contains abundances (read counts) for
# each sequence in all the datasets
sequence_abundance_table <- read.delim(tidy_table_path, sep="\t", header=T)

# make a dataframe containing the sequences and their sequence number (which matches the blast query column)
# and add a column listing their length
sequences <- sequence_abundance_table %>% 
  group_by(sequence, sequence_number) %>% 
  summarize(.groups = "drop") %>%
  mutate(sequence_length = str_length(sequence)) %>% 
  ungroup()

# read in the output file from blastn
blast_df_all <- read.delim(blast_output_path, sep="\t", header=T)

# rename some of the oddly-named columns
blast_df_all <- blast_df_all %>% rename(query = qaccver,
                                        subject = saccver,
                                        percent_identity = pident,
                                        alignment_length = length)


# read in the sample metadata file
metadata_df <- read_excel(sample_metadata_file)

# take the highest scoring blast hit for each query
blast_df <- blast_df_all %>% group_by(query) %>% arrange(-bitscore, .by_group = T) %>% filter(row_number() == 1)

# merge the sequence info and the blast output
blast_df <- left_join(blast_df, sequences, by=c("query" = "sequence_number"))

# what fraction of the query was present in the alignment?
# adding gaps to sequence_length allows gaps in the query relative to the refseq
# TODO: make sure gaps not too high a fraction of the alignment?
blast_df <- blast_df %>% mutate(percent_query_aligned = 100 * alignment_length / (sequence_length + gaps))

# identify sequences that aligned best to internal positive control sequences
#
# these are positive control amplicons that should come up in every sample
# these sequences are marked by their name beginning with INTERNAL_CTRL_
blast_df <- blast_df %>% 
  mutate(internal_control_hit = str_detect(subject, "^INTERNAL_CTRL_"))
         


# This function performs all filtering on datasets and putative assignments.  
#
# The filters for accepting or rejecting an assignment are:
#
#
# Sequence based filters (independent of which datasets contain this sequence): 
#
# 1. Are the sequences similar enough, based on a blastn alignment to the closest
#    reference sequence, and does this alignment cover enough of the sequence length?
#
# Dataset based filters:
# 
# 2. Does the dataset have enough reads in total?
#
# 3. Does the dataset have enough reads mapping to the positive control amplicon?
#
# Sequence x dataset filters (unique for each sequence/dataset combo):
#
# 4. Are their enough reads supporting this sequence in this dataset?
# 
#
filter_assignments <- function(blast_df, sequence_abundance_table, sequences) {
  
  # --------------------------------------------------------------------------
  # logic for saying that one of the observed sequences is close enough 
  # to one of the reference sequences to assign it to that reference sequence
  # --------------------------------------------------------------------------
  
  # thresholds for reference sequence (intended target) assignment
  # these cutoffs are modifiable and are set this way based on empirical 
  # observations of alignments using real datasets
  
  refseq_min_percent_identity <- 98
  refseq_min_percent_aligned <- 99
  
  # thresholds for internal positive control assignment
  internal_ctrl_min_percent_identity <- 90
  internal_ctrl_min_percent_aligned <- 98
  
  # set thresholds based on whether internal control or not
  # since different cutoffs allowed for those two categories.
  blast_df <- blast_df %>% 
    mutate(min_percent_identity = if_else(internal_control_hit, 
                                          internal_ctrl_min_percent_identity, 
                                          refseq_min_percent_identity),
           
           min_percent_aligned = if_else(internal_control_hit, 
                                         internal_ctrl_min_percent_aligned, 
                                         refseq_min_percent_aligned),
           
           # assigned_to_target will be TRUE if this particular sequence
           # aligns to the relevant reference sequence above the thresholds
           assigned_to_target = if_else( (percent_identity > min_percent_identity & 
                                            percent_query_aligned > min_percent_aligned),
                                         TRUE, 
                                         FALSE))
  
  # plot hits: percent identity vs. fraction aligned
  ggplot(blast_df) +
    geom_point(aes(x=percent_query_aligned, 
                   y=percent_identity, 
                   fill=assigned_to_target),
               shape=21, alpha=0.8, size=2) +
    theme_bw() + 
    facet_grid(~internal_control_hit) + 
    xlab("% of query aligned to reference sequence") +
    ylab("% identity between query and reference sequence") +
    # xlim(c(98,100)) +
    # ylim(c(95,100)) +
    scale_fill_manual(values=c("red", "darkslateblue")) +
    {}
  
  # -------------------------------------------
  # separate assigned and unassigned sequences
  # -------------------------------------------
  # assigned sequences, with information about which reference sequence they 
  # matched to, and how closely they did so
  assigned_sequences <- blast_df %>% filter(assigned_to_target)
  
  assigned_queries <- assigned_sequences %>% select(query)
  
  # unassigned sequences are those that weren't assigned to any of the internal control
  # or expected reference sequences
  unassigned_sequences <- filter(sequences, !(sequence_number %in% assigned_sequences$query))
  
  # write out the sequences in fasta format
  # these will be dealt with separately by downstream processes
  writeFasta(unassigned_sequences, "unassigned_sequences.fasta")
  

  # --------------------------------------------
  # merge assignments and abundance information
  # --------------------------------------------

  assigned_sequences_pertinent_columns <- assigned_sequences %>% 
    select(query, subject, percent_identity, gaps, percent_query_aligned, internal_control_hit)
  
  # inner join here discards rows corresponding to unassigned sequences
  # shorter name (sat) for convenience
  sat <- inner_join(sequence_abundance_table, assigned_sequences_pertinent_columns, by=c("sequence_number" = "query"))

  # ----------------------------------------------------------------------------
  # filter out datasets that don't contain a sufficient # of reads aligning to the 
  # internal control positive control sequence
  # ----------------------------------------------------------------------------
  
  # calculate total # of read per dataset and
  #           total # of reads assigned to all internal control targets combined
  # the conditional if_else statement in this mutate allows us to sum the
  # number of reads in each dataset that map to one of the internal control
  # targets
  dataset_summaries <- sat %>% 
    mutate(iab = if_else(internal_control_hit, abundance, 0L)) %>%
    group_by(dataset) %>% 
    summarize(total_dataset_reads = sum(abundance),
              internal_control_reads = sum(iab))  
  
  
  # ------------------------------------------------------
  # filter criterion: minimum total # of reads in dataset
  # ------------------------------------------------------
  
  # the lower cutoff for minimum total # of reads:
  # 0.2 x the mean # of reads in all datasets 
  # this could be modified
  minimum_total_dataset_reads <- mean(dataset_summaries$total_dataset_reads) * 0.2
  
  # ---------------------------------------------------------------------
  # filter criterion: minimum # of reads mapping to internal pos. control
  # ---------------------------------------------------------------------
  
  # the lower cutoff for minimum total # of reads:
  # 0.2 x the mean # of reads in all datasets 
  # this could be modified
  minimum_internal_control_reads <- mean(dataset_summaries$internal_control_reads) * 0.2
  
  dataset_summaries <- dataset_summaries %>% 
    mutate(enough_reads = if_else(total_dataset_reads > minimum_total_dataset_reads, T, F))
  
  dataset_summaries <- dataset_summaries %>% 
    mutate(enough_internal_control_reads = if_else(internal_control_reads > minimum_internal_control_reads, T, F))
  
  
  # TODO: calculate max bin height for histograms (for plotting cutoff)
  # see: https://community.rstudio.com/t/geom-histogram-max-bin-height/10026
  
  # histogram of total reads per dataset 
  ggplot(dataset_summaries) + 
    geom_histogram(aes(x=total_dataset_reads, 
                       fill = enough_reads),
                   bins=60, 
                   color="black")  +
    scale_fill_manual(values = c("grey60", "darkslateblue")) +
    theme_bw() +
    xlab("Total reads in dataset") +
    ylab("Datasets") +
    annotate("segment", x = minimum_total_dataset_reads, xend = minimum_total_dataset_reads, y = 0, yend = 4,
             colour = "red", linetype="dotted", alpha=0.9) +
    
    {}
  
  # histogram of internal_control_abundance values
  ggplot(dataset_summaries) + 
    geom_histogram(aes(x=internal_control_reads, 
                       fill = enough_internal_control_reads),
                   bins=120, 
                   color="black")  +
    scale_fill_manual(values = c("grey60", "darkslateblue")) +
    theme_bw() +
    xlab("Reads mapping to internal positive control (actin) in dataset") +
    ylab("Datasets") +
    annotate("segment", x = minimum_total_dataset_reads, xend = minimum_total_dataset_reads, y = 0, yend = 4,
             colour = "red", linetype="dotted", alpha=0.9) +
    {}
  
  
  # filter out assignments without sufficient reads in datasets
  filtered_datasets <- dataset_summaries %>% 
    filter(enough_reads & enough_internal_control_reads)
  
  # create a subset of the main sequence abundance table based on dataset filtering 
  filtered_sat <- sat %>% filter(dataset %in% filtered_datasets$dataset)
  
  # how many reads needed to make particular assignments?
  # based on an empirical examination of the real data,
  # very little evidence of low-read-count false positive abundances
  # so set a simple minimum # of supporting reads as a cutoff
  #
  # this is a threshold that could be modified or reconsidered
  #
  minimum_supporting_reads_needed <- 20
  
  filtered_sat <- filtered_sat %>% 
    filter(abundance > minimum_supporting_reads_needed)
  
  # look at only non-internal positive control assignments
  assigned_sat <- filtered_sat %>% filter(!internal_control_hit)
  
  # parse out the species of the assignment based on the subject name from the reference sequence
  # which must be of the form NNNNN_NNNNNN (i.e. genus_species)
  # TODO: check this worked for all ref seqs
  assigned_sat <- assigned_sat %>% mutate(species = str_match(subject, "([^_]*_[^_]*)_" )[,2])
  
  # to return multiple objects, must put them in a list:
  # see: https://stackoverflow.com/questions/8936099/returning-multiple-objects-in-an-r-function
  returnList <- list("assigned_sat" = assigned_sat)
  
  return (returnList)
  
}

# perform filtering
filtered_list <- filter_assignments (blast_df, sequence_abundance_table, sequences) 

# pull returned objects out of list
filtered_assignments <- filtered_list$assigned_sat

# make a reporting dataset
# should contain these columns:
#
# *) dataset
# *) target species identified
# *) number of supporting reads
# *) the observed sequence
# *) reference sequence on which assignment was based
# *) percent identity of alignment to reference sequence'
# *) percent of observed sequence contained in the alignment

reporting_dataset <- filtered_assignments %>% select(dataset, species, abundance, sequence, subject, percent_identity, percent_query_aligned)

reporting_dataset <- reporting_dataset %>% 
  rename(
    target_species = species,
    number_supporting_reads = abundance,
    reference_sequence = sequence,
    reference_sequence_name = subject,
    percent_alignment_identity = percent_identity,
    percent_sequence_aligned = percent_query_aligned
)

# write data as plain-text
write.table(reporting_dataset, "identified_targets.tsv", quote=F, sep="\t", col.names=T, row.names=F)

# write as excel
# TODO: date in filename.  Unique run identifier?
wb <- createWorkbook("identified_targets.xlsx")
addWorksheet(wb, "refseq_all")
writeData(wb, "refseq_all", reporting_dataset)
# TODO: save to results dir

# collapse all hits for same target species in each dataset for 
#  include as an additional tab 

reporting_dataset_with_diversity_calcs <- reporting_dataset %>%
	group_by(dataset, target_species) %>% 
	mutate(
	       # pi here is not the number pi, but the proportion of reads from the ith refseq
	       # see: https://en.wikipedia.org/wiki/Diversity_index#Shannon_index
	       pi = number_supporting_reads / sum(number_supporting_reads)) %>%
  ungroup()

reporting_dataset_by_spp <- reporting_dataset_with_diversity_calcs %>% 
	group_by(dataset, target_species) %>% 
  summarize(
            total_number_supporting_reads = sum(number_supporting_reads),   
            average_percent_alignment_identity = mean(percent_alignment_identity),
            average_percent_sequence_aligned = mean(percent_sequence_aligned),
            richness = n(),
            # Shannon index: https://en.wikipedia.org/wiki/Diversity_index#Shannon_index
            # note log() by default calculates natural log
            Shannon_index = -sum(pi*log(pi)),
            # here, evenness is Pielou's evenness index
            # see: https://en.wikipedia.org/wiki/Species_evenness
            # if_else because dividing by log(1) = 0 causes error, so 
            # when richness = 1, evenness = 0 
		        evenness = if_else(richness == 1, 0, Shannon_index / log(richness)), .groups="drop") %>%
  select(-Shannon_index)

addWorksheet(wb, "refseq_collapsed")
writeData(wb, "refseq_collapsed", reporting_dataset_by_spp)


# make a version with rows as datasets and columns as targets.
# include all datasets 
all_datasets <- sequence_abundance_table$dataset



saveWorkbook(wb, "identified_targets.xlsx", overwrite = TRUE)

# TODO: all datasets, wide format, warning flag if something wrong!

# TODO: write the sample metadata to this main output file as a tab





# wider format for a table
# reporting_dataset %>% select(-reference_sequence, -subject, -percent_alignment_identity, -percent_sequence_aligned) %>% pivot_wider(names_from = target_species, values_from = number_supporting_reads) 
# reporting_dataset %>% select(-reference_sequence, -subject, -percent_alignment_identity, -percent_sequence_aligned) %>% pivot_wider(names_from = target_species, values_from = number_supporting_reads) 
wide_df <- reporting_dataset %>% select(-percent_alignment_identity, -percent_sequence_aligned) %>% pivot_wider(names_from = target_species, values_from = number_supporting_reads) 
?pivot_wider




# make a data.table
# ------------------------------------------
# data table with filtering, shading, etc.
# ------------------------------------------
# this will make variant frequencies colored like a heatmap 

# 0 -> 1 by 0.1
# TODO: fix breaks
breaks <- seq(0, 1000, 100)
# create red, green, and blue color scales using RGB encoding
blue_colors <- round(seq(255, 125, length.out = length(breaks) + 1), 0) %>%
  {paste0("rgb(" , . , "," , . , ",255)")}
red_colors  <- round(seq(255, 125, length.out = length(breaks) + 1), 0) %>%
  {paste0("rgb(255," , . , "," , . , ")")}
green_colors  <- round(seq(255, 125, length.out = length(breaks) + 1), 0) %>%
  {paste0("rgb(", . , ", 255 , " , . , ")")}


# create the data table object
dt <- DT::datatable(select(wide_df, -reference_sequence), 
                    caption = 'Identified target species.',
                    filter = 'top',
                    options = list(
                      autoWidth = TRUE,
                      pageLength = 250,
                      # fillContainer = F,
                      # scrollY = TRUE,
                      # scrollX = TRUE,
                      lengthMenu = c(10, 20, 50, 100, 200)
                      # columnDefs = list(list(width = '20px', targets = "_all" ))
                    )) %>%
  # the tail here omits first 4 columns from coloring scheme
  formatStyle(tail(names(wide_df), -3), backgroundColor = styleInterval(breaks, green_colors)) 

dt


