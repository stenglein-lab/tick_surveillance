# stenglein-lab/tick_surveillance pipeline parameters

A pipeline to analyze amplicon sequencing data for tick-borne pathogens

## Input options

Define where the pipeline should find input data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fastq_pattern` | fastq files must match this pattern to be used as input to the pipeline. | `string` | *_R[12]_*.fastq* | True |  |
| `fastq_dir` | A directory containing fastq files that will be input to the pipeline. | `string` |  | True |  |
| `metadata` | The path to an excel spreadsheet containing sample metadata | `string` |  | True |  |
| `refseq_dir` | Default location of certain input files | `string` | refseq/ |  |  |
| `targets` | A file containing the target sequences and other information about these sequences. | `string` | ${refseq_dir}/targets.tsv | True |  |

## Primer and adapter trimming options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `primers` | A file containing the names and sequences of primers used to amplify targets. | `string` | ${refseq_dir}/primers.tsv | True |  |
| `post_trim_min_length` | After trimming of adapter and primer sequences, amplicons shorter than this length will be discarded. | `integer` | 100 | True |  |
| `amplicon_primers_max_error_fraction` | This specifies the error tolerance (fraction) used when searching for adapters sequences to trim. This value is passed to the the cutadapt -e parameter. | `number` | 0.2 | True |  |
| `adapters_min_overlap` | Used in trimming of Illumina adapter sequences.  Specifies the minimum length of overlap between adapter sequence and read sequence for trimming to occur.  Passed to cutadapt -O parameter. | `integer` | 10 | True |  |

## Calling of positives and Surveillance Report



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `surveillance_columns` | This file specifies what columns will be included in the surveillance report and optional default values. | `string` | ${refseq_dir}/surveillance_columns.txt | True |  |
| `min_reads_for_positive_surveillance_call` | The number of read pairs assigned to a particular target in order for that target to be called positive. | `integer` | 50 | True |  |
| `max_blast_refseq_evalue` | The maximum BLASTN e-value for an ASV to be initially assigned to a target sequence.  Final assignment will be based on additional criteria. | `number` | 1e-10 | True |  |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  |  |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  |  |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  |  |

## DADA2 options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `dada_outdir` | Directory in which DADA2 output files will be placed. | `string` | ${outdir}/dada2 | True |  |

## Tree-building

Parameters related to the construction of phylogenetic trees that include observed sequences and pre-defined reference sequences.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `make_trees` | A flag to turn on or off this step. | `boolean` | True | True |  |
| `tree_outdir` | A directory in which tree files will be placed. | `string` | ${outdir}/trees | True |  |

## BLASTing of unassigned sequences

Parameters associated with BLASTing of unassigned sequences against the NCBI nt database for the purposes of taxonomic classification of non-target sequences.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `blast_unassigned_sequences` | A flag to turn off this BLAST step | `boolean` |  |  |  |
| `local_nt_database_dir` | The path to a directory containing a local copy of the NCBI nt database (or whatever database you wish to search).  Note that this is not the path of the database, but is the path of the directory containing the database. | `string` |  |  |  |
| `local_nt_database_name` | The name of the local BLAST database to be searched in this step. | `string` | nt |  |  |
| `blast_tax_dir` | The directory containing local copies of the files contained in https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz.  This is optional and if you don't specify it (recommended), the pipeine will download these automatically. | `string` |  |  |  |
| `remote_blast_nt` | Setting this parameter to true will run this BLAST search using NCBI's remote copy of the nt blast database.  This avoids needing to have a local copy of the nt database installed.  Note that this will be very slow. | `boolean` |  |  |  |
| `max_blast_nt_evalue` | The maximum BLAST e-value for classifying unassigned sequences. | `number` | 1e-10 | True |  |
| `blast_perc_identity` | The minimum percent identity for BLAST hits to be considered when classifying unassigned sequences. | `number` | 70 | True |  |
| `blast_qcov_hsp_perc` | The minimum query coverage percentage for BLAST hits to be considered when classifying unassigned sequences. | `number` | 70 | True |  |

## Software dependencies



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `singularity_pull_docker_container` | If you are using singularity and are persistently observing issues downloading Singularity images directly due to timeout or network issues, then you can use the --singularity_pull_docker_container parameter to pull and convert the Docker image instead. | `boolean` |  |  | True |
| `R_tar_dir` | Directory containing R package .tar.gz files | `string` | ${baseDir}/lib/R/ |  | True |
| `python_venv_path` | The path to a directory where the python venv will be created | `string` | python_venv/ |  | True |
| `python_requirements` | A file listing the python packages that will be installed in the venv | `string` | ${baseDir}/lib/requirements.txt |  | True |
| `script_dir` | The path to a directory where scripts for creating python venv exist | `string` | ${baseDir}/scripts/ |  | True |

## Output options

Output related options

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `outdir` | The output directory where results will be saved. You have to use absolute paths to store on Cloud infrastructure. | `string` | results/ | True |  |
| `output_prefix` | Main pipeline output file names will contain this prefix. | `string` | Today's date in YYYY_MM_DD format | True |  |
| `blast_outdir` | Directory in which BLAST output will be placed | `string` | ${outdir}/blast | True |  |
| `trimmed_outdir` | Directory in which adapter, primer, and quality trimmed fastq will be placed | `string` | ${outdir}/trimmed_fastq | True |  |
| `QC_and_summary_stats_dir` | Directory in which QC report files will be placed and other files summarizing pipeline run and output | `string` | ${outdir}/QC_and_summary_stats |  |  |
| `multiQC_reports` | Directory in which multiqc reports will be placed | `string` | ${QC_and_summary_stats_dir}/multiqc |  |  |
| `tracedir` | Directory in which Nextflow logs and reports will be placed. | `string` | ${QC_and_summary_stats_dir}/pipeline_info |  |  |
| `log_outdir` | Directory in which pipeline logs and reports will be placed. | `string` | ${QC_and_summary_stats_dir}/log |  |  |
| `cutadapt_trim_reports` | Directory in which cutadapt trim reports will be placed. | `string` | ${QC_and_summary_stats_dir}/cutadapt_trim_reports |  |  |
| `publish_dir_mode` | The default mode for publishing files in output dirs (hard link) | `string` | link |  |  |

## Unassigned sequence filtering



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `filter_unassigned_seq` | List of taxa of interest | `string` | Borrelia,Borreliella,Babesia,Anaplasma,Ehrlichia |  |  |

## Reporting options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `email` | Email address that will be mailed on pipeline completion | `string` |  |  |  |
| `email_on_fail` | Email address that will be mailed on pipeline completion in case of failure | `string` |  |  |  |
| `plaintext_email` | Non-HTML email only | `boolean` |  |  | True |
| `monochrome_logs` | Don't use color in logs | `boolean` |  |  | True |
| `multiqc_title` |  | `string` |  |  | True |
| `multiqc_config` |  | `string` |  |  | True |
| `multiqc_logo` |  | `string` |  |  | True |
| `max_multiqc_email_size` |  | `string` | 25.MB |  | True |
| `multiqc_methods_description` |  | `string` |  |  | True |

## Parameter validation and display

Parameters related to parameter validation and display

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `validate_params` | Validate parameter values using this schema | `boolean` |  |  | True |
| `schema_ignore_params` | Ignore validation for this list of parameters | `string` |  |  | True |
| `show_hidden_params` | Show parameters marked as hidden in this schema when running pipeline | `boolean` |  |  | True |
| `validationLenientMode` | Validate parameters in lenient mode.  See: https://nextflow-io.github.io/nf-validation/parameters/validation/#variable-type-checking | `boolean` | True |  |  |
