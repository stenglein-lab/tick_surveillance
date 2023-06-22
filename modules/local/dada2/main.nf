/*
  This process runs dada2 on input files
*/
process DADA2 {      
  publishDir "${params.dada_outdir}", mode: 'link'
  tag        "all"

  label 'process_high'
  label 'error_retry'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0"
  } else {
      container "quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0"
  }

  input:
  path(all_trimmed_reads)

  output:
  path("dada_seqtab.txt")                    , emit: seqtab
  path("dada_filtered")                      , emit: dada_filtered_dir
  path "versions.yml"                        , emit: versions                                         

  script:
  """
  # Run dada2 using trimmed fastq as input and create a tabular output of results
  Rscript ${params.script_dir}/run_dada_on_trimmed.R .
  """
}

process TIDY_DADA_OUTPUT {      
  publishDir "${params.dada_outdir}", mode: 'link'
  tag        "all"

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }

  input:
  path(dada_seqtab) 

  output:
  path "observed_sequences.fasta"      , emit: observed_sequences
  path "sequence_abundance_table.tsv"  , emit: sequence_abundance_table
  path "versions.yml"                  , emit: versions                                         

  script:
  """
  # This R script creates an output named observed_sequences.fasta containing all of the
  # unique sequences observed in the amplicon dataset
  # and a sequence_abundance_table.tsv, which lists the abundances of these
  # sequences in each dataset
  Rscript ${params.script_dir}/tidy_dada_output.R $dada_seqtab
  """
}

