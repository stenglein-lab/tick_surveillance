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

# if running from Rscript (called from the pipeline)
if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  r_bindir = args[1]
  tidy_table_path = args[2]
  blast_output_path = args[3]
  sample_metadata_file = args[4]
  targets_tsv_file = args[5]
  output_dir = "./"
} else {
  # if running via RStudio
  r_bindir  =  "."
  tidy_table_path = "../results/sequence_abundance_table.tsv"
  blast_output_path = "../results/observed_sequences.fasta.bn_refseq"
  sample_metadata_file = "../input/AK_metadata.xlsx"
  targets_tsv_file = "../refseq/targets.tsv"
  output_dir = "../results/"
}


# a function to write out sequences in fasta format
# assumes the dataframe has columns named sequence_number and observed_sequence 
# see: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"sequence_number"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"observed_sequence"]))
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

# make sure the dataset names have character types, in case they have integer value names (1, 2, ...)
sequence_abundance_table$dataset <- as.character(sequence_abundance_table$dataset)

# make a dataframe containing the sequences and their sequence number (which matches the blast query column)
# and add a column listing their length
sequences <- sequence_abundance_table %>% 
  group_by(sequence, sequence_number) %>% 
  summarize(.groups = "drop") %>%
  mutate(sequence_length = str_length(sequence)) %>% 
  rename(observed_sequence = sequence) %>% 
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
                                        alignment_length = length,
                                        query_length = qlen,
                                        subject_length = slen)

# take the highest scoring blast hit for each query
blast_df <- blast_df_all %>% group_by(query) %>% arrange(-bitscore, .by_group = T) %>% filter(row_number() == 1)

# TODO: Deal with the case when a sequence hits multiple refseqs equally well 
# This would be the case when a sequence is equally close to 
# two of the reference sequences.  E.g. a Borrelia burgdorferi equally similar to B31 and N40
blast_df_equal_hits <- blast_df_all %>% group_by(query) %>% mutate(max_bitscore = max(bitscore)) %>% filter(bitscore == max_bitscore) %>% summarize(n=n()) %>% filter(n>1)

# merge the sequence info and the blast output
blast_df <- left_join(blast_df, sequences, by=c("query" = "sequence_number")) 

# 
# what fraction of the query was present in the alignment?
# 
# Scenarios for insertions or deletions in the alignment:
# 
# 1) If observed sequences has no indels relative to the reference
#   - alignment_length == query length
#   - gaps = 0
#
# 2) If observed sequences has a deletion relative to the reference
#   - alignment_length == (query_length + gaps)
#   - gaps = # of bases in deletion
#
# 3) If observed sequences has an insertion relative to the reference
#   - alignment_length == query length
#   - gaps = # of bases in insertion
#
# Note that any indels relative to reference sequence have the effect of decreasing 
# the % identity of the alignment in blast output in addition to influencing 
# gaps and alignment length values as described above.  So, any gaps in the alignment
# will produce a lower-scoring alignment less likely to be above cutoff.
#  

blast_df <- blast_df %>% mutate(percent_query_aligned = 100 * alignment_length / query_length,
                                percent_of_alignment_gaps = 100 * gaps / alignment_length)

# ----------
# metadata 
# ----------

# read in the sample metadata file
metadata_df <- read_excel(sample_metadata_file)

# keep track of original metadata names
original_metadata_names <- names(metadata_df)

# fix the metadata names using tidyverse make.names function
names(metadata_df) <- make.names(names(metadata_df),unique = TRUE)

# Make sure the Index column is character data, in case Indexes have integer values (1, 2, etc.)
metadata_df$Index <- as.character(metadata_df$Index)

# assign datasets to batches
if ( ! "batch" %in% colnames(metadata_df)) {
  # if no batch column specified in the metadata,
  # then consider all of the datasets as belonging to a single batch
  message ('INFO: There is no "batch" column in the metadata file so will consider all datasets to belong to a single batch')
  metadata_df <- metadata_df %>% mutate(batch = "1")
}

# ---------------------
# target information
# ---------------------
targets_df <- read.delim(targets_tsv_file, sep="\t", header=T)
targets_df$internal_control <- as.logical(targets_df$internal_control)

# rename target sequences to be clear that they include primer sequences, 
# which are stripped off of the actual observed sequences 
# since primer-derived bases are unreliable
targets_df <- targets_df %>% rename(target_sequence_incl_primers = sequence)

# convert empty strings in reporting column to NA values
targets_df$reporting <- na_if(targets_df$reporting, "")

# join target info with blast info
# TODO: is this join necessary?   
blast_df <- left_join(blast_df, targets_df, by=c("subject" = "ref_sequence_name"))

# --------------------------------------------------------------------------
# logic for saying that one of the observed sequences is close enough 
# to one of the reference sequences to assign it to that reference sequence
# --------------------------------------------------------------------------
assign_blast_hits_to_refseqs <- function() {
  
  # <<- is the global assignment operator 
  # (it modifies the global blast_df variable instead of creating a new variable 
  # with the same name in the function scope)
  blast_df <<- blast_df %>%  mutate(

    # assigned_to_target will be TRUE if this particular sequence
    # aligns to the relevant reference sequence above the thresholds
    # the min_percent_identity and min_percent_aligned come from the targets_tsv_file
    # so can be assigned separately for each target
    assigned_to_target = if_else( (percent_identity > min_percent_identity & 
                                     percent_query_aligned > min_percent_aligned &
                                     percent_of_alignment_gaps < max_percent_gaps),
                                  TRUE, 
                                  FALSE))
}  


# Call function to see if sequences are close enough to a target to be assigned to that target
assign_blast_hits_to_refseqs()

# plot hits: percent identity vs. fraction aligned
assigned_plot <- ggplot(blast_df) +
  geom_point(aes(x=percent_query_aligned, 
                 y=percent_identity, 
                 fill=assigned_to_target),
             shape=21, alpha=0.8, size=2) +
  scale_fill_manual(values=c("red", "darkslateblue")) +
  theme_bw() + 
  facet_wrap(~subject) +
  xlab("% of query aligned to reference sequence") +
  ylab("% identity between query and reference sequence") 

# render the plot
# TODO: output this as a plot/report ?
assigned_plot



# --------------------------------------------------------------------------------------
# separate assigned and unassigned sequences and write out unassigned sequences in fasta
# --------------------------------------------------------------------------------------
# assigned sequences
assigned_sequences <- blast_df %>% filter(assigned_to_target)

# unassigned sequences are those that weren't assigned to any target reference sequence
unassigned_sequences <- filter(sequences, !(sequence_number %in% assigned_sequences$query))

# write out the unassigned sequences in fasta format
# these will be dealt with separately by downstream processes
writeFasta(unassigned_sequences, paste0(output_dir, "unassigned_sequences.fasta"))

# -----------------------
# consolidate dataframes
# -----------------------

# get rid of 0 counts
sparse_sat <- filter(sequence_abundance_table, abundance > 0) 

# keep track of metadata rows
metadata_key <- metadata_df %>% select(Index, batch)

# join sparse SAT with blast_df 
dataset_df <- full_join(metadata_key, sparse_sat, by=c("Index" = "dataset"))
dataset_df <- left_join(dataset_df, blast_df, by=c("sequence_number" = "query")) %>% filter(assigned_to_target == T)

# ---------------------------------------------------------------------
# QC criterion: minimum # of reads mapping to internal pos. control
# this corresponds to the Acceptable DNA column in the reporting data
# ---------------------------------------------------------------------

# calculate the # of tick actin mapping reads in real tick (non-control) datasets
non_control_dataset_batch_averages <- dataset_df %>%
  filter(internal_control) %>%
  group_by(batch, Index) %>%
  mutate(per_sample_internal_ctrl_reads = log10(sum(abundance))) %>%
  ungroup() %>%
  group_by(batch) %>%
  summarize(mean_batch_internal_ctrl_reads = mean(per_sample_internal_ctrl_reads),
            sd_batch_internal_ctrl_reads = sd(per_sample_internal_ctrl_reads))

# join in info about batch averages 
dataset_df <- left_join(dataset_df, non_control_dataset_batch_averages, by="batch")

# cutoffs
dataset_df <- dataset_df %>%
  mutate(minimum_internal_control_log_reads =
           mean_batch_internal_ctrl_reads - (3 * sd_batch_internal_ctrl_reads),
         # right now, set minimum # of reads for real targets at 10...
         minimum_non_control_reads = 10)


# -------------------------------------------------------
# collapse unique sequences to level of targeted species
# -------------------------------------------------------

# collapse all hits for same target species in each dataset 
dataset_by_spp <- dataset_df %>% 
  group_by(Index, species) %>%
  mutate(
    abundance = sum(abundance),   
    percent_identity = mean(percent_identity),
    percent_query_aligned = mean(percent_query_aligned),
    richness = n(),
  ) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  # get rid of unneeded columns
  select(Index, batch, species, abundance, percent_identity, 
         percent_query_aligned, richness, internal_control, 
         minimum_internal_control_log_reads, minimum_non_control_reads)

tick_reads_p <- dataset_df %>% 
  filter(internal_control) %>% 
  group_by(batch, Index, species) %>%
  summarize(abundance = sum(abundance), above_cutoff = log10(abundance) > minimum_internal_control_log_reads) %>%
  ggplot() +
  geom_violin(aes(x=batch, y=abundance), color="slateblue", fill="NA", size=0.5) +
  geom_jitter(aes(x=batch, y=abundance, fill=above_cutoff),
              shape=21, size=3, stroke=0.15, width=0.2) +
  scale_fill_manual(values=c("red", "darkslateblue")) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Reads mapping to tick actin") +
  xlab("Batch") +
  ggtitle("Number of tick-mapping reads in individual tick or control datasets\n")

tick_reads_p


# ----------
# make calls
# ----------
dataset_df_calls <- dataset_by_spp %>% 
  mutate(pos_neg_call = case_when((internal_control & (log10(abundance) >= minimum_internal_control_log_reads)) ~ "Positive", 
                                  (!internal_control & (abundance >= minimum_non_control_reads)) ~ "Positive", 
                                  TRUE ~ "Negative"),
         .groups="drop") %>%
  select(Index, species, abundance, pos_neg_call)



# ------------------------------------
# fill in missing samples and targets
# ------------------------------------

# dataset_df_calls is a sparse dataset at this point
# it only contains rows corresponding to datasets that had any reads matching one of the targets 
# so we need to fill in the missing datasets: those datasets with no assigned targets
# and those species that were not observed in any dataset

# fill out the dataframe with any species listed in the targets file but 
# not present in the data (species for which there are no reads in any dataset)
df_species <- targets_df %>% group_by(species)  %>% summarize()

df_species_observed  <- filter(dataset_df_calls, pos_neg_call == "Positive") %>% group_by(species) %>% summarize() %>% filter(!is.na(species))

df_species_not_observed <- as.data.frame(df_species$species[!(df_species$species %in% df_species_observed$species)])
colnames(df_species_not_observed) <- c("species")

# make zeroed out rows for species not observed
num_extra_spp <- nrow(df_species_not_observed)

new=1
if (new) {
# create data frame rows for species that didn't appear in any dataset
# what if we don't have any calls?  
# situation possible with just a few datasets
if (nrow(dataset_df_calls) == 0){ 
  extra_spp_rows_one <- data.frame(Index = character(),
                                   species = character(),
                                   abundance = integer(),
                                   pos_neg_call = character(),
                                   stringsAsFactors=FALSE)
} else {
  # get the first row of the dataset_df_calls data frame as a place holder for the unobserved species
  extra_spp_rows_one <- head(dataset_df_calls, 1)
}

# replicate this row as many times as their are extra species
extra_spp_rows <- slice_sample(extra_spp_rows_one, n=num_extra_spp, replace=T)
} else {
  extra_spp_rows <- head(dataset_df_calls, num_extra_spp)
}

# fill in a negative call for each unobserved species
extra_spp_rows$species <- df_species_not_observed$species
extra_spp_rows$abundance <- rep(0L, num_extra_spp)
extra_spp_rows$pos_neg_call <- rep("Negative", num_extra_spp)

# add on these extra rows corresponding to 0 counts for species not observed in data
dataset_df_calls <- rbind(dataset_df_calls, extra_spp_rows)

# add in samples that had no data (no positive calls).  do this by left_join() to metadata
dataset_df_reporting <- left_join(metadata_key, dataset_df_calls) %>% select(-c(batch, abundance))

# pivot to a wider format for reporting
reporting_dataset <- dataset_df_reporting %>% 
  pivot_wider(names_from=species, 
              values_from=c(pos_neg_call),
              names_sep = "_")

# pivot to a wider format for reporting
# add in samples that had no data (no positive calls).  do this by left_join() to metadata
dataset_df_reporting_abundances <- left_join(metadata_key, dataset_df_calls) %>% select(-c(batch, pos_neg_call))

reporting_dataset_abundances <- dataset_df_reporting_abundances %>% 
  pivot_wider(names_from=species, 
              values_from=c(abundance),
              names_sep = "_")

# ----------------------------
# create a surveillance table
# ----------------------------

# output date in format specified by CDC surveillance folks
# R/tidyverse happier with these column names without spaces or characters like "-" or "(" or ")"
surv_cols_ok_names <- 
             c("Index", 
               "ID", "Submission_No", "Vial_ID", "Pathogen_Testing_ID", 
               "CSID", "No_Samples_in_Testing_Pool", "Extraction_SOP_Number", 
               "Extraction_SOP_Version", "Acceptable_DNA", 
               "Borrelia_sp", "Borrelia_burgdorferi_sensu_stricto", 
               "Borrelia_mayonii", "Borrelia_miyamotoi", 
               "Borrelia_Other_species_name", 
               "Anaplasma_phagocytophilum_strain_not_differentiated", 
               "Babesia_microti", "Ehrlichia_muris_eauclairensis", 
               "Powassan_virus_lineage_I", "DeerTick_virus_Powassan_virus_lineage_II", 
               "Attempted_to_Sequence", "Successfully_Sequenced", "In_house_Tick_Molecular_ID", 
               "In_house_Blood_Meal_ID", "In_house_Tick_Host_Molecular_ID", 
               "In_house_Vertebrate_Molecular_ID", "Notes")

# how big should the DF be?
num_surv_cols <- length(surv_cols_ok_names)
num_surv_rows <- nrow(reporting_dataset)  

# create an empty DF
surv_df <- data.frame(mat = matrix(ncol = num_surv_cols, nrow = num_surv_rows), stringsAsFactors = F)

# name the columns
colnames(surv_df) <- surv_cols_ok_names

# -----------------------------------------------
# populate the columns in the surveillance table
# -----------------------------------------------

# --------
# metadata
# --------

# Illumina Index
surv_df$Index <- reporting_dataset$Index

# Pull in metadata from the metadata dataframe (from excel input)
populate_surveillance_metadata <- function(surv_df, meta_df, column_names) {
  
  for (column_name in column_names) {
    if (!(column_name %in% colnames(meta_df))) {
      message (paste0("WARNING: column ", column_name, " not present in metadata file: ", sample_metadata_file, ".  Will not populate column in surveillance report table."))
      next
    }
    surv_df[column_name] <- meta_df[column_name]
  }
  return (surv_df)
}

metadata_columns <-  c("ID", "Submission_No", "Vial_ID", "Pathogen_Testing_ID", 
                       "CSID", "No_Samples_in_Testing_Pool", "Extraction_SOP_Number", 
                       "Extraction_SOP_Version")

surv_df <- populate_surveillance_metadata(surv_df, metadata_df, metadata_columns)
               
# ------------------------
# Positive/Negative Calls
# ------------------------

# acceptable DNA (reads mapping to tick actin)
surv_df$Acceptable_DNA <- replace_na(reporting_dataset$Ixodes_scapularis, "Negative")
# swith Pos/Neg -> True/False for Acceptable DNA column
surv_df$Acceptable_DNA <- recode(surv_df$Acceptable_DNA, Positive = "TRUE", Negative="FALSE")

# Borrelia burgdorferi
surv_df$Borrelia_burgdorferi_sensu_stricto <- replace_na(reporting_dataset$Borrelia_burgdorferi, "Negative")

# Borrelia mayonii
surv_df$Borrelia_mayonii <- replace_na(reporting_dataset$Borrelia_mayonii, "Negative")

# Borrelia miyamotoi
surv_df$Borrelia_miyamotoi <- replace_na(reporting_dataset$Borrelia_miyamotoi, "Negative")

# Anaplasma phagocytophilum
surv_df$Anaplasma_phagocytophilum_strain_not_differentiated <- replace_na(reporting_dataset$Anaplasma_phagocytophilum, "Negative")

# Babesia microti
surv_df$Babesia_microti <- replace_na(reporting_dataset$Babesia_microti, "Negative")

# EMLA
# Ehrlichia_muris_eauclairensis", 
surv_df$Ehrlichia_muris_eauclairensis <- replace_na(reporting_dataset$Ehrlichia_EMLA, "Negative")

# Deal with other Borrelia spp
borrelia_sp_df <- reporting_dataset %>% 
  select(Borrelia_andersoni, Borrelia_bissettii) %>%
  mutate(
    
    # Positive for other Borrelia spp?
    Borrelia_sp = case_when(
    # if either of these spp are positive then other Borrelia_sp is also pos.
    Borrelia_andersoni == "Positive" ~ "Positive",
    Borrelia_bissettii  == "Positive" ~ "Positive",
    TRUE ~ "Negative" ),
    
    # specify other Borrelia spp names
    # this is awkward.  Consider improving 
    Borrelia_Other_species_name = case_when(
    Borrelia_andersoni == "Positive" & Borrelia_bissettii == "Positive" ~ "Borrelia andersonii,Borrelia bissettiae",

    Borrelia_andersoni == "Positive" ~ "Borrelia andersonii",
    Borrelia_bissettii  == "Positive" ~ "Borrelia bissettiae",
    TRUE ~ NA_character_)
  ) 

# assign to columns in main surveillance df
surv_df$Borrelia_sp <- replace_na(borrelia_sp_df$Borrelia_sp, "Negative")
surv_df$Borrelia_Other_species_name <- borrelia_sp_df$Borrelia_Other_species_name

# not tested columns for POWV
surv_df$Powassan_virus_lineage_I <- rep("Not Tested", num_surv_rows)
surv_df$DeerTick_virus_Powassan_virus_lineage_II <- rep("Not Tested", num_surv_rows)


# make a version of the surveillance dataframe with abundance info instead of +/- calls
surv_df_abundances <- surv_df

# acceptable DNA (reads mapping to tick actin)
surv_df_abundances$Acceptable_DNA <- replace_na(reporting_dataset_abundances$Ixodes_scapularis, 0)

# Borrelia burgdorferi
surv_df_abundances$Borrelia_burgdorferi_sensu_stricto <- replace_na(reporting_dataset_abundances$Borrelia_burgdorferi, 0)

# Borrelia mayonii
surv_df_abundances$Borrelia_mayonii <- replace_na(reporting_dataset_abundances$Borrelia_mayonii, 0)

# Borrelia miyamotoi
surv_df_abundances$Borrelia_miyamotoi <- replace_na(reporting_dataset_abundances$Borrelia_miyamotoi, 0)

# Anaplasma phagocytophilum
surv_df_abundances$Anaplasma_phagocytophilum_strain_not_differentiated <- replace_na(reporting_dataset_abundances$Anaplasma_phagocytophilum, 0)

# Babesia microti
surv_df_abundances$Babesia_microti <- replace_na(reporting_dataset_abundances$Babesia_microti, 0)

# EMLA
# Ehrlichia_muris_eauclairensis", 
surv_df_abundances$Ehrlichia_muris_eauclairensis <- replace_na(reporting_dataset_abundances$Ehrlichia_EMLA, 0)

# Deal with other Borrelia spp
borrelia_sp_abundances_df <- reporting_dataset_abundances %>% 
  select(Borrelia_andersoni, Borrelia_bissettii) %>%
  mutate( Borrelia_sp = sum(Borrelia_andersoni, Borrelia_bissettii, na.rm=T))

surv_df_abundances$Borrelia_sp <- replace_na(borrelia_sp_abundances_df$Borrelia_sp, 0)


# -------------------
# Write out results
# -------------------

# write data as plain-text
write.table(dataset_df, paste0(output_dir, "identified_targets.tsv"), quote=F, sep="\t", col.names=T, row.names=F)

# write as excel
# TODO: date in filename?  Unique run identifier?

wb <- createWorkbook(paste0(output_dir, "identified_targets.xlsx"))
modifyBaseFont(wb, fontSize = 11, fontColour = "black", fontName = "Helvetica")

# a generic style for all cells
all_cell_style <- createStyle( 
  border = "TopBottomLeftRight",
  borderColour = getOption("openxlsx.borderColour", "grey"),
  borderStyle = getOption("openxlsx.borderStyle", "thin"),
  halign = "left",
  numFmt = "0.00",
  wrapText = F
)

# integer number format
integer_num_style <- createStyle( 
  numFmt = "0"
)

# column headers
col_header_style <- createStyle( 
  fontName = "Helvetica",
  fontSize = 11,
  textDecoration = "bold",
  wrapText = F
)


# this function will do some basic styling for worksheets
style_worksheet <- function (wb, sheetname, df) {
  
  # style all cells
  addStyle(wb=wb, sheet=sheetname,
           style=all_cell_style,
           cols = 1:ncol(df)+1,
           rows = 1:nrow(df)+1,
           gridExpand = T
  )
  
  # style column headers
  addStyle(wb=wb, sheet=sheetname,
           style=col_header_style,
           cols = 1:ncol(df)+1,
           rows = 1:1,
           gridExpand = T,
           stack = T
  )
  
  # set column widths to auto
  setColWidths(wb, sheet = sheetname, cols = 1:ncol(df), widths = "auto")
}

# populate the workbook with data

addWorksheet(wb, "surveillance")
writeData(wb, "surveillance", surv_df)
style_worksheet(wb, "surveillance", surv_df)

addWorksheet(wb, "surveillance_counts")
writeData(wb, "surveillance_counts", surv_df_abundances)
style_worksheet(wb, "surveillance_counts", surv_df_abundances)
addStyle(wb, "surveillance_counts", integer_num_style, rows=1:nrow(surv_df_abundances)+1, cols=9:ncol(surv_df_abundances)+1, gridExpand =T, stack = T)

addWorksheet(wb, "data_by_species")
writeData(wb, "data_by_species", dataset_by_spp)
style_worksheet(wb, "data_by_species", dataset_by_spp)

addWorksheet(wb, "all_data")
writeData(wb, "all_data", dataset_df)
style_worksheet(wb, "all_data", dataset_df)

addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata_df)
style_worksheet(wb, "metadata", metadata_df)

addWorksheet(wb, "targets")
writeData(wb, "targets", targets_df)
style_worksheet(wb, "targets", targets_df)


# add conditional formatting for surveillance worksheet

# light green fill
light_green_fill <- createStyle( 
  bgFill = "#CCFF99"
)

# red fill
red_fill <- createStyle( 
  bgFill = "#FF9999"
)

# Acceptable DNA column 
conditionalFormatting(wb=wb, sheet="surveillance", 
                      "colourScale",
                      cols = 10:10,
                      rows = 1:nrow(surv_df)+1,
                      style = light_green_fill,
                      rule = "TRUE",
                      type = "contains"
)

# Acceptable DNA column 
conditionalFormatting(wb=wb, sheet="surveillance", 
                      "colourScale",
                      cols = 10:10,
                      rows = 1:nrow(surv_df)+1,
                      style = red_fill,
                      rule = "FALSE",
                      type = "contains"
)

# Pos/Neg calls - red color for positive calls
conditionalFormatting(wb=wb, sheet="surveillance", 
                      "colourScale",
                      cols = 11:18,
                      rows = 1:nrow(surv_df)+1,
                      style = red_fill,
                      rule = "Positive",
                      type = "contains"

                      )

# Read counts
conditionalFormatting(wb=wb, sheet="surveillance_counts", 
                      "colourScale",
                      cols = 10:18,
                      rows = 1:nrow(surv_df_abundances)+1,
                      style = light_green_fill,
                      rule = ">0",
                      type = "expression"
)
# Pos/Neg calls - red color for positive calls
conditionalFormatting(wb=wb, sheet="surveillance", 
                      "colourScale",
                      cols = 11:18,
                      rows = 1:nrow(surv_df)+1,
                      style = red_fill,
                      rule = "Positive",
                      type = "contains"
)


# write out the workbook
saveWorkbook(wb, paste0(output_dir, "sequencing_report.xlsx"), overwrite = TRUE)


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

