# Tick-borne pathogen surveillance bioinformatics pipeline                    

This repository contains a bioinformatics pipeline for the analysis of amplicon sequencing datasets of tick-associated microbes.  This was developed by the Stenglein lab at Colorado State University in collaboration with researchers in CDC's Division of Vector-Borne Diseases.

This pipeline is [described in this paper](https://pubmed.ncbi.nlm.nih.gov/37247570/).

## Running the pipeline

See the [dependencies section](#dependencies) below for information about the main dependencies required for running this pipeline(including nextflow and singularity).

### Running from github

The simplest way to run the pipeline is directly from github, like this:

```
nextflow run stenglein-lab/tick_surveillance -resume --metadata /path/to/metadata_xls --fastq_dir /path/to/fastq/directory -profile singularity
```

You must specify two required inputs to the pipeline: the path to a metadata excel spreadsheet and the path to a directory containing input fastq.  See [this section](#inputs) for more information on required inputs.

### Running test datasets

The pipeline includes a handful of small datasets (<= 1000 read) that are derived from real known positive (or known negative) datasets.  These are included in the [test directory](./test/) of the repository.  These datasets serve as positive and negative controls and allow you to test that the pipeline is working as expected.  To use these test datasets, run with the test profile, for instance:

```
nextflow run stenglein-lab/tick_surveillance -profile singularity,test
```

Or to run with conda:
```
nextflow run stenglein-lab/tick_surveillance -profile conda,test
```

##### Test outputs:

The results of the test run will be placed in a `test/results` sub-directory.

#### Running a specific version from github:

To run a specific version of the pipeline, use the -r option, for example:

```
nextflow run stenglein-lab/tick_surveillance -profile singularity,test -r v1.0.7
```

#### Making sure that you are running the latest version when running from github.

Nextflow will cache the pipeline in `$HOME/.nextflow/assets/` and continue to use the cached version, even if the pipeline has newer versions on github.  To remove the locally cached version, which will force nextflow to download and cache the latest version, run:

```
nextflow drop stenglein-lab/tick_surveillance
# now run with latest version
nextflow run stenglein-lab/tick_surveillance -profile singularity,test
```

Alternatively, you can just delete the cached pipeline directory:
```
rm -rf ~/.nextflow/assets/stenglein-lab/tick_surveillance/
```

Running from github is [described in more detail here](https://www.nextflow.io/docs/latest/sharing.html).  



#### Running in different environments

You will want to use a profile that matches your computing environment.  So, for instance, if running on an SGE HPC environment, you'd run something like:

```
nextflow run stenglein-lab/tick_surveillance -resume --metadata /path/to/metadata_xls --fastq_dir /path/to/fastq/directory -profile singularity,sge 
```

### Running by cloning the pipeline's repository

It is also possible to download the pipeline code to a directory of your choosing.  This can be useful if, for instance, you want to modify or debug the code.  You can do this by cloning the repository (or a fork of the repository):

```
git clone https://github.com/stenglein-lab/tick_surveillance.git
cd tick_surveillance
nextflow run main.nf -resume --metadata /path/to/metadata_xls --fastq_dir /path/to/fastq/directory -profile singularity
```

## Inputs

The pipeline requires two inputs:

1. [A metadata file](#metadata-file).

2. [Sequence datasets in fastq format](#input-fastq).  

### Metadata file

A metadata file in Microsoft Excel format must be input to the pipeline.  

1. The spreadsheet should contain a single tab
2. The single tab should contain a table that defines metadata for samples being analyzed.
3. The first row of the table should contain column names.
4. Subsequent rows should contain metadata with one row per dataset. 
5. One of the columns must be named Index (case sensitive), and the values in this column must match fastq file names.  In other words, the values in this column should match the sample names specified in the Illumina sample sheet.  

An example of a working metadata file [can be found here](./test/test_metadata.xlsx)

### Input fastq

Input sequence data is assumed to be Illumina paired-end data in separate read1 and read2 files.  Files can be compressed or not but it would be preferred to leave them as compressed files to save disk space.

The expected names of the fastq files are defined by the parameter `fastq_pattern`, whose default value is defined in nextflow.config as `*_R[12]_*.fastq*`.  This pattern can be overridden on the nextflow command line using the `--fastq_pattern` parameter.

The location of the fastq files is specified by the `fastq_dir` parameter, whose default value is `./input/fastq` (relative to your present working directory from which you are running the nextflow command).

It is expected that sample IDs are not repeated in the Illumina sample sheet.

## Output

The main outputs of the pipeline are:

1. [A surveillance report](#surveillance-report)
2. [QC reports](#qc-reports)
3. Information about [observed sequences that were not assigned to any reference sequences](#unassigned-sequences).

#### Output directory

Output files are placed in a `results` directory (or `test/results` when running with -profile test).  The output directory can be specified using the `--outdir` parameter (e.g. `nextflow run stenglein-lab/tick_surveillance ... --outdir a_results_directory`

#### Output file name prefixes 

The main pipeline output file names will be prefixed by a value that is by default the date the pipeline is run (e.g. `2023_04_06_sequencing_report.xlsx`).  This filename prefix can be changed using the --output_prefix parameter.  For instance, running `nextflow run stenglein-lab/tick_surveillance ... --output_prefix my_new_run` will create a file named `my_new_run_sequencing_report.xlsx`)

### Surveillance Report

The pipeline outputs a surveillance report in Microsoft Excel format.

#### Surveillance columns

The columns in this report are defined in [this file](./refseq/surveillance_columns.txt).  It is possible to add or remove columns from this report by adding or removing them from this file.  

This is a 2-column tab-delimited file.  The first column contains the names of the columns that will form the surveillance report table.  The second column contains optional default text for thiis column.

To add, remove columns from the surveillance report table, add or remove lines from this file.  Columns can also be reordered by reordering lines in this file.

#### Data in surveillance report

The values in the surveillance report come from two possible sources:

1. **Metadata**.  

2. **Read counts from the sequence data.**

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
| reporting_columns      | A semicolon-separated list of column names in the reporting table.  Reads that are assigned to this reference sequence will be assigned to these columns in the surveillance table. |
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
nextflow run stenglein-lab/tick_surveillance -profile singularity --primers /path/to/primers.tsv
```

If primer sequences are not entered in the correct orientation, trimming will not work and targets amplified by these sequences will not be detected by the pipeline.  The solution in this case will most likely just be to swap the F/R orientation of the primers in this file.  

**Correct primer orientation:**  The forward primer (primer_f) should appear at the beginning of Illumina read 1 and be in the same orientation as R1.  The reverse primer (primer_r) should appear at the beginning of read 2 and be in the same orientation as read 2.  In other words, relative to the PCR product as a whole, the primers should point towards each other.

## Dependencies

This pipeline has two main dependencies: nextflow and singularity.  These programs must be installed on your computer to run this pipeline.

### Nextflow

To run the pipeline you will need to be working on a computer that has nextflow installed. Installation instructions can be [found here](https://www.nextflow.io/docs/latest/getstarted.html#installation).  To test if you have nextflow installed, run:

```
nextflow -version
```

This pipeline requires nextflow version > 22.10

### Singularity 

The pipeline uses singularity containers to run programs like cutadapt, BLAST, and R.  To use these containers you must be running the pipeline on a computer that has [singularity](https://sylabs.io/singularity) [installed](https://sylabs.io/guides/latest/admin-guide/installation.html).  To test if singularity is installed, run:

```
singularity --version
```

Singularity containers will be automatically downloaded and stored in a directory named `singularity_cacheDir` in your home directory.  They will only be downloaded once.

### Conda

It is possible to run this pipeline using an all-in-one [conda](https://docs.conda.io/en/latest/) environment, defined [in this file](./conda/tick_conda_environment.yaml).  But it is strongly recommended to use singularity instead of conda.  

### R libraries

Some of the pipeline code is implemented in [R scripts](./scripts/).  Some of these scripts require R packages like [openxlsx](https://www.rdocumentation.org/packages/openxlsx/versions/4.2.5.2), for writing output in Excel format.  These packages are installed locally, on top of a Rocker tidyverse singularity image.  This occurs in nextflow process `setup_R_dependencies`, which invokes [this script](./scripts/install_R_packages.R).

### Python dependencies

Some of the pipeline code is implemented in [Python scripts](./scripts/).  In particular, the tree-building scripts.  These python scripts require various python modules.  This is handled by creating a python virtual environment (venv), which happens in nextflow process `setup_python_venv`.  This venv is then activated from a basic python singularity image (for instance in process `create_fasta_for_trees`).

## BLASTing of unassigned sequences

It is possible that amplicon sequencing will generate sequences that are off-target or not closely related enough to be assigned to one of the predefined reference sequences.  The pipeline can BLAST these "unassigned" sequences against the NCBI nt database to try to figure out what they are.  

Enabling BLASTing of unassigned sequences is controlled by the `blast_unassigned_sequences` parameter.  

There are two ways to run the BLAST:

1. **Remote option:** By using BLASTn with the -remote option.  This is very slow but doesn't require a local copy of the BLAST database.  To run this way, set the `remote_blast_nt` parameter to true.  This sends the sequences to a remote NCBI server for BLASTing.
2. **Local option:** Alternatively, if you have a local copy of the nt BLAST database installed, you can specify its location using the `local_nt_database_dir` parameter, which should be the path to a directory containing a local nt database.  The name of this database is expected to be "nt", but this name can be changed by overriding the `local_nt_database_name` parameter.

These options can be configured on the command line, for example:

```
nextflow run stenglein-lab/tick_surveillance -profile test,singularity --blast_unassigned_sequences --remote_blast_nt true
```

Or in a nextflow config file, for instance:

```
  // profile for local BLAST of unassigned sequences
  local_blast {
    params.blast_unassigned_sequences = true
    params.local_nt_database_dir ="/home/NCBI_databases/nt/"
    params.remote_blast_nt = false
  }
```

#### BLAST Taxonomy 

BLAST can provide taxonomic information about database hits.  The pipeline downloads these files automatically as part of each run.


## Parameter information


## Primer and adapter trimming options                                                                                            
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `primers` | A file containing the names and sequences of primers used to amplify targets. | `string` | refseq/primers.tsv | True
| `post_trim_min_length` | After trimming of adapter and primer sequences, amplicons shorter than this length will be discarded. |
| `amplicon_primers_max_error_fraction` | This specifies the error tolerance (fraction) used when searching for adapters sequences
| `adapters_min_overlap` | Used in trimming of Illumina adapter sequences.  Specifies the minimum length of overlap between adapte
                                                                                                                                  
## Calling of positives and Surveillance Report                                                                                   
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `surveillance_columns` | This file specifies what columns will be included in the surveillance report and optional default value
| `min_reads_for_positive_surveillance_call` | The number of read pairs assigned to a particular target in order for that target t
| `max_blast_refseq_evalue` | The maximum BLASTN e-value for an ASV to be initially assigned to a target sequence.  Final assignme
                                                                                                                                  
## Input/output options                                                                                                           
                                                                                                                                  
Define where the pipeline should find input data and save output data.                                                            
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `fastq_pattern` | fastq files must match this pattern to be used as input to the pipeline. | `string` | *_R[12]_*.fastq* | True 
| `fastq_dir` | A directory containing fastq files that will be input to the pipeline. | `string` | input/fastq/ | True |  |      
| `outdir` | The output directory where results will be saved. You have to use absolute paths to store on Cloud infrastructure. | 
| `targets` | A file containing the target sequences and other information about these sequences. | `string` | refseq//targets.tsv
| `tracedir` | Directory in which Nextflow logs and reports will be placed. | `string` | ${params.outdir}/pipeline_info |  |  |   
| `log_outdir` | Directory in which pipeline logs and reports will be placed. | `string` | ${params.outdir}/log |  |  |           
| `output_prefix` | Main pipeline output file names will contain this prefix. | `string` | Today's date in YYYY_MM_DD format | Tru
                                                                                                                                  
## Max job request options                                                                                                        
                                                                                                                                  
Set the top limit for requested resources for any single job.                                                                     
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set 
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to 
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set 
                                                                                                                                  
## Generic options                                                                                                                
                                                                                                                                  
Less common options for the pipeline, typically set in a config file.                                                             
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `help` | Display help text. | `boolean` |  |  |  |                                                                              
                                                                                                                                  
## DADA2 options                                                                                                                  
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `dada_outdir` | Directory in which DADA2 output files will be placed. | `string` | ${params.outdir}/dada2 | True |  |           
                                                                                                                                  
## BLASTing of unassigned sequences                                                                                               
                                                                                                                                  
Parameters associated with BLASTing of unassigned sequences against the NCBI nt database for the purposes of taxonomic classificat
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `blast_unassigned_sequences` | A flag to turn off or on this BLAST / assignment step. | `boolean` |  | True |  |                
| `max_blast_nt_evalue` | The maximum BLAST e-value. | `number` | 1e-10 | True |  |                                               
| `blast_perc_identity` | The minimum percent identity for BLAST hits to be considered. | `number` | 70 | True |  |               
| `blast_qcov_hsp_perc` | The minimum query coverage percentage for BLAST hits to be considered. | `number` | 70 | True |  |      
| `local_nt_database_dir` | The path to a directory containing a local copy of the NCBI nt database (or whatever database you wish
| `local_nt_database_name` | The name of the local BLAST database to be searched in this step. | `string` | nt | True |  |        
| `remote_blast_nt` | Setting this parameter to true will run this BLAST search using NCBI's remote copy of the nt blast database.
| `blast_tax_dir` | The directory containing local copies of the files contained in https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.ta
                                                                                                                                  
## Tree-building                                                                                                                  
                                                                                                                                  
Parameters related to the construction of phylogenetic trees that include observed sequences and pre-defined reference sequences. 
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `make_trees` | A flag to turn on or off this step. | `boolean` | True | True |  |                                               
| `tree_outdir` | A directory in which tree files will be placed. | `string` | ${params.outdir}/trees | True |  |                 
                                                                                                                                  
## Other parameters                                                                                                               
                                                                                                                                  
| Parameter | Description | Type | Default | Required | Hidden |                                                                  
|-----------|-----------|-----------|-----------|-----------|-----------|                                                         
| `python_venv_path` |  | `string` | python_venv/ |  | True |                                                                     
| `script_dir` |  | `string` | scripts/ |  |  |                                                                                   
| `initial_fastqc_dir` |  | `string` | ${params.outdir}/initial_fastqc/ |  |  |                                                   
| `post_trim_fastqc_dir` |  | `string` | ${params.outdir}/post_trim_fastqc/ |  |  |                                               
| `trimmed_outdir` |  | `string` | ${params.outdir}/trimmed_fastq |  |  |                                                         
| `blast_outdir` |  | `string` | ${params.outdir}/blast |  |  |                                                                   
| `refseq_dir` |  | `string` | /home/mdstengl/2023_1_issues/refseq/ |  |  |                                                       
| `singularity_pull_docker_container` |  | `boolean` |  |  |  |                                                                   
| `metadata` |  | `string` | None |  |  |                                                                                         
| `publish_dir_mode` |  | `string` | link |  |  |                                                                                 
| `multiqc_title` |  | `string` |  |  |  |                                                                                        
| `multiqc_config` |  | `string` |  |  |  |                                                                                       
| `multiqc_logo` |  | `string` |  |  |  |                                                                                         
| `max_multiqc_email_size` |  | `string` | 25.MB |  |  |                                                                          
| `multiqc_methods_description` |  | `string` |  |  |  |                                                                          
| `email` |  | `string` | mark.stenglein@colostate.edu |  |  |                                                                    
| `email_on_fail` |  | `string` | mark.stenglein@colostate.edu |  |  |                                                            
| `plaintext_email` |  | `string` |  |  |  |                                                                                      
| `monochrome_logs` |  | `string` |  |  |  |                                                                                      
| `hook_url` |  | `string` |  |  |  |                                                                                             
| `version` |  | `string` |  |  |  |                                                                                              
| `validate_params` |  | `string` |  |  |  |                                                                                      
| `show_hidden_params` |  | `string` |  |  |  |                                                                                   
| `python_requirements` |  | `string` | /home/mdstengl/2023_8_dsl2_dev/lib/requirements.txt |  |  |                               
| `QC_and_summary_stats` |  | `string` | results/QC_and_summary_stats |  |  |                                                     
| `cutadapt_trim_reports` |  | `string` | results/QC_and_summary_stats/cutadapt_trim_reports |  |  |                              
| `multiQC_reports` |  | `string` | results/QC_and_summary_stats/multiqc |  |  |                                                  
| `schema_ignore_params` |  | `string` |  |  |  |                                                                                 
| `R_tar_dir` |  | `string` | /home/mdstengl/2023_8_dsl2_dev/lib/R/ |  |  |                                                       
| `R_install_pkg_script` |  | `string` | /home/mdstengl/2023_8_dsl2_dev/bin/install_R_packages.R |  |  |                          
| `validate_metadata_script` |  | `string` | /home/mdstengl/2023_8_dsl2_dev/bin/validate_metadata.R |  |  |                       
                                                                                                                                  



