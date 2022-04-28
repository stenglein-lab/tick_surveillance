# Tick-borne pathogen surveillance bioinformatics pipeline                    

This repository contains a bioinformatics pipeline for the analysis of amplicon sequencing datasets of tick-associated microbes

For more information, please see [the documentation](./documentation)


## Input to the pipeline

The pipeline requires two main inputs:

1. A metadata file.

2. Sequence data in fastq format.

## Surveillance Report

The pipeline outputs a surveillance report in Microsoft Excel format.  

#### Surveillance columns

The columns in this report are defined in [this file](./refseq/surveillance_columns.txt).  It is possible to add or remove columns from this report by adding or removing them from this file.  

This is a 2-column tab-delimited file.  The first column contains the names of the columns that will form the surveillance report table.  The second column contains optional default text for thiis column.

To add, remove columns from the surveillance report table, add or remove lines from this file.  Columns can also be reordered by reordering lines in this file.

#### Data in surveillance report

The values in the surveillance report come from two possible sources:

1. **Metadata**.  

2. **Read counts from the sequence data.**


## To add a new reference sequence (a new target)

Reference sequence (aka targets) are defined in the [targets.tsv](./refseq/targets.tsv) file.  This tab-delimited file contains the following columns: 

| Column                 | Description |
| -----------            | ----------- |
| ref_sequence_name      | The reference sequecne name |
| species                | The species for this reference sequence.  This will be reported but is not used to populate the surveillance table (reporting_columns is used for that). |
| reporting_columns      | A comma-separate list of column names in the reporting table.  Reads that are assigned to this reference sequence will be assigned to these columns in the surveillance table. |
| min_percent_identity   | The minimum percent identity of the alignment between an observed sequence and this reference sequence to be assigned as a positive hit |
| min_percent_aligned    | The minimum percent of the observed sequence that must align to this reference sequence to be assigned as a positive hit |
| max_percent_gaps       | The maximum percent gap characters in alignments of observed sequences and this reference sequence to be assigned as a positive hit|
| internal_control       | True if this corresponds to an internal control target, such as tick actin |
| sequence               | The expected reference sequence, including primers |

Note that the names of the reporting_columns specified in this file must match exactly the names of columns defined in the surveillance columns definition file.

#### To add a new reference sequence (a new target)

To add a new sequence to the targets.tsv file, you will need to edit this file.  It is a plain-text [tab-delimited file](https://en.wikipedia.org/wiki/Tab-separated_values) that can be edited in google sheets or similar software.  Add a new row for the new target sequence.  From google sheets, download the file in tab-separated value format and transfer it to the computer where you will be running this pipeline.  

The default location of the targets.tsv file can be overriden by specifying the --targets option on the nextflow command line.  For instance:

```
nextflow run main.nf -profile singularity,CDC_hpc --targets /path/to/targets.tsv
```

## To add new primer pairs

Primers are defined in the [primers.tsv](./refseq/primers.tsv) file.  Primer sequences defined in this file are used for two purposes:

1. To identify read pairs that have expected forward and reverse primers at their ends in the expected F/R orientation.  Only read pairs with a pair of expected primers in the expected orientation will be kept for further analysis.
2. Primer sequences will be trimmed off of observed sequences, since primer-derived sequences do not reliably reflect the template sequence (in case of primer-template mismatches).

To add a new primer pair to the pipeline, you will need to edit this file.  It is a plain-text [tab-delimited file](https://en.wikipedia.org/wiki/Tab-separated_values) that can be edited in google sheets or similar software.  Add a new row for the new primer pair.  From google sheets, download the file in tab-separated value format and transfer it to the computer where you will be running this pipeline.  

The default location of the primers.tsv file can be overwritten by specifying the --primers option on the nextflow command line.  For instance:

```
nextflow run main.nf -profile singularity,CDC_hpc --primers /path/to/primers.tsv
```

If primer sequences are not entered in the correct orientation, trimming will not work and targets amplified by these sequences will not be detected by the pipeline.  The solution in this case will most likely just be to swap the F/R orientation of the primers in this file. 
