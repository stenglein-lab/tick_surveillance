/*
  This process prepends main output files with a prefix (by default the date, can be overridden)
*/

process PREPEND_OUTPUT_FILENAMES {
  publishDir "${params.outdir}", mode: 'link'

  label 'process_low'

  input:
  path(output_file) 

  output:
  path ("${params.output_prefix}${output_file}"), emit: prepended_file

  script:
  """
  cp $output_file ${params.output_prefix}${output_file}
  """
}

