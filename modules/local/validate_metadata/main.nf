/*
  This process validates that the metadata file is in the approprate
  format and contains the appropriate information.  Logic in the R script.
*/
process VALIDATE_METADATA {
  tag        "${metadata}"
  publishDir "${params.outdir}", mode: "copy"


  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }

  input:
  path (validate_script) 
  path (metadata) 
  path (sample_ids) 

  output:
  path(metadata)       , emit: validated_metadata
  path(sample_ids)     , emit: sample_ids
  path "versions.yml"  , emit: versions                                         


  script:
  """
    Rscript ${validate_script} $metadata $sample_ids
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
  """
}

