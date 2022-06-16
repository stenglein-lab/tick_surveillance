# Tick-borne pathogen surveillance bioinformatics pipeline                    

This repository contains a bioinformatics pipeline for the analysis of amplicon sequencing datasets of tick-associated microbes.  This was developed by the Stenglein lab at Colorado State University in collaboration with researchers in CDC's Division of Vector-Borne Diseases.


## Input 

The pipeline requires two inputs:

1. A [metadata file](#metadata-file).

2. [Sequence datasets in fastq format](#input-fastq).  

### Metadata file

A metadata file in Microsoft Excel format must be input to the pipeline.  

1. The spreadsheet should contain a single tab
2. The single tab should contain a table that defines metadata for samples being analyzed.
3. The first row of the table should contain column names.
4. Subsequent rows should contain metadata with one row per dataset. 
5. One of the columns must be named Index (case sensitive), and the values in this column must match fastq file names.  In other words, the values in this column should match the sample names specified in the Illumina sample sheet.

TODO: Provide a metadata template.

### Input fastq

Input sequence data is assumed to be Illumina paired-end data in separate read1 and read2 files.  Files can be compressed or not but it would be preferred to leave them as compressed files to save disk space.

The expected names of the fastq files are defined by the parameter `fastq_pattern`, whose default value is defined in nextflow.config as `*_R[12]_*.fastq*`.  This pattern can be overridden on the nextflow command line using the `--fastq_pattern` paramter.

The location of the fastq files is specified by the `fastq_dir` parameter, whose default value is `${baseDir}/input/fastq`.

## Output

The main outputs of the pipeline are:

1. [A surveillance report](#surveillance-report)
2. [QC reports](#qc-reports)
3. Information about [observed sequences that were not assigned to any reference sequences](#unassigned-sequences).

### Surveillance Report

The pipeline outputs a surveillance report in Microsoft Excel format in the results directory (or value of parameter `outdir`).   

#### Surveillance columns

The columns in this report are defined in [this file](./refseq/surveillance_columns.txt).  It is possible to add or remove columns from this report by adding or removing them from this file.  

This is a 2-column tab-delimited file.  The first column contains the names of the columns that will form the surveillance report table.  The second column contains optional default text for thiis column.

To add, remove columns from the surveillance report table, add or remove lines from this file.  Columns can also be reordered by reordering lines in this file.

#### Data in surveillance report

The values in the surveillance report come from two possible sources:

1. **Metadata**.  

2. **Read counts from the sequence data.**

TODO: flesh out this section.

### QC reports

Two QC reports in HTML format from initial input sequence data (`...initial_qc_report.html`) and from sequence data after trimming low quality and adapter bases (`...post_trim_qc_report.html`) are output.

### Unassigned sequences

Observed sequences that were not assigned to any reference sequence are used as queries to a BLASTN search against the NCBI nt database.  The output of this blast search is summarized in a `...non_reference_sequence_assignments.xlsx` output file: an MS Excel spreadsheet.

## Reference sequences

### To add a new reference sequence (a new target)

Reference sequence (aka targets) are defined in the [targets.tsv](./refseq/targets.tsv) file.  This tab-delimited file contains the following columns: 

| Column                 | Description |
| -----------            | ----------- |
| ref_sequence_name      | The reference sequecne name |
| species                | The species for this reference sequence.  This will be reported but is not used to populate the surveillance table (reporting_columns is used for that). |
| primer_name            | The name of the primers expected to amplify this target.  Provided for reference only.
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
nextflow run main.nf -profile singularity --targets /path/to/targets.tsv
```

## To add new primer pairs

Primers are defined in the [primers.tsv](./refseq/primers.tsv) file.  Primer sequences defined in this file are used for two purposes:

1. To identify read pairs that have expected forward and reverse primers at their ends in the expected F/R orientation.  Only read pairs with a pair of expected primers in the expected orientation will be kept for further analysis.
2. Primer sequences will be trimmed off of observed sequences, since primer-derived sequences do not reliably reflect the template sequence (in case of primer-template mismatches).

To add a new primer pair to the pipeline, you will need to edit this file.  It is a plain-text [tab-delimited file](https://en.wikipedia.org/wiki/Tab-separated_values) that can be edited in google sheets or similar software.  Add a new row for the new primer pair.  From google sheets, download the file in tab-separated value format and transfer it to the computer where you will be running this pipeline.  

The default location of the primers.tsv file can be overwritten by specifying the --primers option on the nextflow command line.  For instance:

```
nextflow run main.nf -profile singularity --primers /path/to/primers.tsv
```

If primer sequences are not entered in the correct orientation, trimming will not work and targets amplified by these sequences will not be detected by the pipeline.  The solution in this case will most likely just be to swap the F/R orientation of the primers in this file. 
