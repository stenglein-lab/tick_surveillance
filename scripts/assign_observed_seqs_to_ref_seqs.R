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

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#

# if running from Rscript
if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  # TODO: check CLAs
  r_bindir = args[1]
  tidy_table_path = args[2]
  blast_output_path = args[3]
  sample_metadata_file = args[4]
  targets_csv_file = args[5]
  output_dir = "./"
} else {
  # if running via RStudio
  r_bindir  =  "."
  tidy_table_path = "../results/sequence_abundance_table.tsv"
  blast_output_path = "../results/observed_sequences.fasta.bn_refseq"
  # sample_metadata_file = "../input/sample_metadata.xlsx"
  sample_metadata_file = "../input/Group2_metadata.xlsx"
  targets_csv_file = "../refseq/targets.csv"
  output_dir = "../results/"
}


# a function to write out sequences in fasta format
# assumes the dataframe has columns named sequence_number and sequence 
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


# ------------------
# read input files
# -----------------

# ----------------------------
# sequence abundance table
# ----------------------------
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

# -------------
# blast input
# -------------

# read in the output file from blastn
blast_df_all <- read.delim(blast_output_path, sep="\t", header=T)

# rename some of the oddly-named columns
blast_df_all <- blast_df_all %>% rename(query = qaccver,
                                        subject = saccver,
                                        percent_identity = pident,
                                        alignment_length = length)

# take the highest scoring blast hit for each query
# TODO: what if a sequence hits multiple refseqs equally well?  
blast_df <- blast_df_all %>% group_by(query) %>% arrange(-bitscore, .by_group = T) %>% filter(row_number() == 1)

# merge the sequence info and the blast output
blast_df <- left_join(blast_df, sequences, by=c("query" = "sequence_number"))

# what fraction of the query was present in the alignment?
# adding gaps to sequence_length allows gaps in the query relative to the refseq
# TODO: make sure gaps not too high a fraction of the alignment?
blast_df <- blast_df %>% mutate(percent_query_aligned = 100 * alignment_length / (sequence_length + gaps))

# ----------
# metadata 
# ----------

# read in the sample metadata file
metadata_df <- read_excel(sample_metadata_file)

# fix the metadata names using tidyverse make.names function
names(metadata_df) <- make.names(names(metadata_df),unique = TRUE)

# identify control samples : what type are they?
metadata_df <- metadata_df %>% mutate(
  control_type = case_when(
    Pathogen.testing.Sample.No. == "EC"  ~ "EC",
    Pathogen.testing.Sample.No. == "CTC" ~ "CTC",
    Pathogen.testing.Sample.No. == "KF"  ~ "KF",
    Pathogen.testing.Sample.No. == "NTC" ~ "NTC",
    Pathogen.testing.Sample.No. == "PC"  ~ "PC",
    TRUE                                 ~ NA_character_ ))

# assign datasets to batches
if ( ! "batch" %in% colnames(metadata_df)) {
  # if no batch column specified in the metadata,
  # assume batch corresponds to the characters of the Index prior to the "-"
  
  batch_pattern <- "(\\S+)-"
  
  # check that we can identify such a pattern
  if (!all(str_detect(metadata_df$Index, batch_pattern))){
    stop ('There is no "batch" column in the metadata file and the "Index" column does not contain batch information in the expected format. Exiting.')
  }
    
  # pull out the batch info from the Index value
  metadata_df <- metadata_df %>% mutate(batch = str_match(metadata_df$Index, batch_pattern)[,2])   
}

# ---------------------
# target information
# ---------------------
targets_df <- read.delim(targets_csv_file, sep=",", header=T)
targets_df$internal_control <- as.logical(targets_df$internal_control)
targets_df$reported <- as.logical(targets_df$reported)

# join target info onto blast_df info
blast_df <- left_join(blast_df, targets_df, by=c("subject" = "ref_sequence_name"))

#
#  Are the sequences similar enough, based on a blastn alignment to the closest
#    reference sequence, and does this alignment cover enough of the sequence length?
#
assign_blast_hits_to_refseqs <- function() {
  
  # --------------------------------------------------------------------------
  # logic for saying that one of the observed sequences is close enough 
  # to one of the reference sequences to assign it to that reference sequence
  # --------------------------------------------------------------------------
  
  # <<- is the global assignment operator 
  # (it modifies the global blast_df variable instead of creating a new variable 
  # with the same name in the function scope)
  blast_df <<- blast_df %>%  mutate(
    # assigned_to_target will be TRUE if this particular sequence
    # aligns to the relevant reference sequence above the thresholds
    # the min_percent_identity and min_percent_aligned come from the targets_csv_file
    # so can be assigned separately for each target
    
    assigned_to_target = if_else( (percent_identity > min_percent_identity & 
                                     percent_query_aligned > min_percent_aligned),
                                  TRUE, 
                                  FALSE))
  
  # plot hits: percent identity vs. fraction aligned
  # TODO: output this as a plot/report
  assigned_plot <- ggplot(blast_df) +
    geom_point(aes(x=percent_query_aligned, 
                   y=percent_identity, 
                   fill=assigned_to_target),
               shape=21, alpha=0.8, size=2) +
    theme_bw() + 
    facet_grid(~internal_control) + 
    xlab("% of query aligned to reference sequence") +
    ylab("% identity between query and reference sequence") +
    scale_fill_manual(values=c("red", "darkslateblue")) +
    {}
  
  # to return multiple objects, must put them in a list:
  # see: https://stackoverflow.com/questions/8936099/returning-multiple-objects-in-an-r-function
  returnList <- list("blast_df" = blast_df, "assigned_plot" = assigned_plot)
}  


# Pull out close but not quite close enough sequences

assign_returns <- assign_blast_hits_to_refseqs()

# plot assigned vs. not assigned
assign_returns$assigned_plot




# -------------------------------------------
# separate assigned and unassigned sequences
# -------------------------------------------
# assigned sequences, with information about which reference sequence they 
# matched to, and how closely they did so
assigned_sequences <- blast_df %>% filter(assigned_to_target)

# unassigned sequences are those that weren't assigned to any of the internal control
# or expected reference sequences
unassigned_sequences <- filter(sequences, !(sequence_number %in% assigned_sequences$query))

# write out the sequences in fasta format
# these will be dealt with separately by downstream processes
writeFasta(unassigned_sequences, paste0(output_dir, "unassigned_sequences.fasta"))

# -----------------------
# consolidate dataframes
# -----------------------

# join metadata w/ sequence abundance table

# get rid of 0 counts
sparse_sat <- filter(sequence_abundance_table, abundance > 0) 

# join sparse SAT with blast_df 
# sparse_sat_with_blast <- left_join(sparse_sat, blast_df, by=c("sequence_number" = "query"))  %>% filter(!is.na(subject))

# inner_join duplicates metadata rows for datasets with >1 target
dataset_df <- inner_join(metadata_df, sparse_sat, by=c("Index" = "dataset")) %>% select(-sequence)

dataset_df <- left_join(dataset_df, blast_df, by=c("sequence_number" = "query")) %>% filter(assigned_to_target == T)



# ----------------------------------------------------------------------------
# filter out datasets that don't contain a sufficient # of reads aligning to the 
# internal control positive control sequence
# ----------------------------------------------------------------------------

# ---------------------------------------------------------------------
# filter criterion: minimum # of reads mapping to internal pos. control
# ---------------------------------------------------------------------

# first, determine how many reads map to the internal control (tick actin) in
# the EC (extraction control) datasets
# do this separately for each batch of datasets (each plate unless specified otherwise)
extraction_control_dataset_summaries <- filter(dataset_df, control_type == "EC" & internal_control) 


# calculate average (mean) tick actin-mapping reads in each batch
#
extraction_control_batch_averages <- extraction_control_dataset_summaries %>%
  group_by(batch) %>% 
  # per Andrias's suggestion: drop max and min EC values for calculating 
  # average of control (tick) reads
  filter(abundance < max(abundance) & abundance > min(abundance)) %>%
  summarize(mean_EC_internal_ctrl_reads = mean(abundance))
            

# this will fill in 0 values for batches without any EC tick actin mapping reads
batches <- metadata_df %>% group_by(batch) %>% summarize()
extraction_control_batch_averages <- left_join(batches, extraction_control_batch_averages, by="batch") %>% 
  mutate(mean_EC_internal_ctrl_reads = replace_na(mean_EC_internal_ctrl_reads, 0))

# use the number of internal control reads in the extraction control datasets
# to set the threshold for the minimum # of internal control reads in the 
# real sample datasets
# 
# If this number is zero, use an alternative minimum cutoff, 
# based on the # tick reads in the real sample datasets

#
# The # of tick actin mapping reads in datasets appears to be
# lognormally distributed
#
# test if distribution is actually lognormal
dataset_test_lognormal <- dataset_df %>%
  filter(is.na(control_type) & internal_control) %>%
  group_by(batch, Index) %>%
  summarize(per_sample_internal_ctrl_reads = log10(sum(abundance))) %>%
  # summarize(per_sample_internal_ctrl_reads = (sum(abundance))) %>%
  ungroup() %>%
  filter(batch == "D")

# test if distribution is actually lognormal
shapiro.test(dataset_test_lognormal$per_sample_internal_ctrl_reads)


# calculate the # of tick actin mapping reads in real tick (non-control) datasets
non_control_dataset_batch_averages <- dataset_df %>%
  filter(is.na(control_type) & internal_control) %>%
  group_by(batch, Index) %>%
  mutate(per_sample_internal_ctrl_reads = log10(sum(abundance))) %>%
  ungroup() %>%
  group_by(batch) %>%
  summarize(mean_batch_internal_ctrl_reads = mean(per_sample_internal_ctrl_reads),
            sd_batch_internal_ctrl_reads = sd(per_sample_internal_ctrl_reads))

# join in info about batch averages 
dataset_df <- left_join(dataset_df, extraction_control_batch_averages, by="batch")
dataset_df <- left_join(dataset_df, non_control_dataset_batch_averages, by="batch")

# Plot # of tick-mapping reads in individual datasets
tick_reads_p <- dataset_df %>% 
  filter(internal_control) %>% 
  group_by(batch, Index, control_type, species) %>%
  summarize(abundance = sum(abundance)) %>%
  ggplot() +
  geom_violin(aes(x=batch, y=abundance), fill="NA", size=0.5) +
  geom_jitter(aes(x=batch, y=abundance, fill=control_type), shape=21, size=3, stroke=0.15, width=0.2) +
  scale_y_log10() +
  scale_fill_manual(values=c("blue", "red", "green", "yellow", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Reads mapping to tick actin") +
  xlab("Plate") +
  ggtitle("Number of tick-mapping reads in individual tick or control datasets\n")

tick_reads_p

# Are there enough tick-mapping reads to call tick samples positive?
# if there are reads in the EC (extraction control) datasets,
# use the mean value of these (excluding max and min)
#
# Otherwise, if no EC tick actin reads,  use three standard deviations below the mean value 
# for the tick datasets
#
# TODO: will have to finalize the rules here.
#
dataset_df <- dataset_df %>% 
  mutate(minimum_internal_control_reads = 
           if_else(
             (mean_EC_internal_ctrl_reads > 0),
             mean_EC_internal_ctrl_reads,
             mean_batch_internal_ctrl_reads - (3 * sd_batch_internal_ctrl_reads)),
         # right now, set minimum # of reads for real targets at 50...
         minimum_non_control_reads = 50)

tick_df <- filter(dataset_df, internal_control & control_type == "EC")


# -------------------------------------------------------
# collapse unique sequences to level of targeted species
# -------------------------------------------------------
# collapse all hits for same target species in each dataset 
dataset_by_spp <- dataset_df %>% group_by(Index, species) %>%
  mutate(
    abundance = sum(abundance),   
    percent_identity = mean(percent_identity),
    percent_query_aligned = mean(percent_query_aligned),
    richness = n(),
  ) %>%
  filter(row_number() == 1)

# ----------
# make calls
# ----------
dataset_df_calls <- dataset_by_spp %>% 
  mutate(pos_neg_call = if_else(internal_control, (abundance > minimum_internal_control_reads), (abundance > minimum_non_control_reads)))

dataset_df_summary <- dataset_df_calls %>% 
  select(-sequence_number, -subject, -mismatch, -gaps, -qstart, -qend, -sstart, -send, 
         -evalue, -bitscore, -sequence_length, -alignment_length,
         -min_percent_identity, -min_percent_aligned, 
         -assigned_to_target, -mean_EC_internal_ctrl_reads, -sequence,
         -sd_batch_internal_ctrl_reads,
         -mean_batch_internal_ctrl_reads, -minimum_internal_control_reads, 
         -percent_query_aligned,
         -minimum_non_control_reads, -richness, -reported, -internal_control)

reporting_dataset <- dataset_df_summary %>% 
  pivot_wider(names_from=species, 
              values_from=c(pos_neg_call, abundance, percent_identity), 
              names_sep = " ")


# write data as plain-text
write.table(reporting_dataset, "identified_targets.tsv", quote=F, sep="\t", col.names=T, row.names=F)

# write as excel
# TODO: date in filename.  Unique run identifier?

wb <- createWorkbook("identified_targets.xlsx")

addWorksheet(wb, "summary")
writeData(wb, "summary", reporting_dataset)

addWorksheet(wb, "all_data_table")
writeData(wb, "all_data_table", dataset_df)

saveWorkbook(wb, "identified_targets.xlsx", overwrite = TRUE)

# TODO: all datasets, wide format, warning flag if something wrong!

# TODO: write the sample metadata to this main output file as a tab


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
# dt <- DT::datatable(select(wide_df, -reference_sequence), 
#                     caption = 'Identified target species.',
#                     filter = 'top',
#                     options = list(
#                       autoWidth = TRUE,
#                       pageLength = 250,
#                       # fillContainer = F,
#                       # scrollY = TRUE,
#                       # scrollX = TRUE,
#                       lengthMenu = c(10, 20, 50, 100, 200)
#                       # columnDefs = list(list(width = '20px', targets = "_all" ))
#                     )) %>%
#   # the tail here omits first 4 columns from coloring scheme
#   formatStyle(tail(names(wide_df), -3), backgroundColor = styleInterval(breaks, green_colors)) 
# 
# dt
