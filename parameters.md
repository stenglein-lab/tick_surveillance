

# stenglein-lab/tick_surveillance pipeline parameters                                                                             
                                                                                                                                  
A pipeline to analyze amplicon sequencing data for tick-borne pathogens                                                           
                                                                                                                                  
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
                                                                                                                                  


