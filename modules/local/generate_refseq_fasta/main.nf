workflow GENERATE_REFSEQ_FASTA {

  take: 
  targets

  main:
  generate_refseq_fastas(targets)
  combine_refseq_fastas(generate_refseq_fastas.out.fasta.collect())
  
  emit:
  fasta = combine_refseq_fastas.out.fasta
  // no specialized tool versions to emit
}

/*
  generate one fasta for each reference sequence
*/
process generate_refseq_fastas {
  label 'process_low'
  tag "$targets.refseq_name"

  input:
  val targets 

  output:
  path ("*fasta"), emit: fasta

  // makes fasta formatted records for each target
  script:
  """
  printf ">%s\n%s\n" $targets.ref_sequence_name $targets.sequence > ${targets.ref_sequence_name}.fasta
  """
}

/*
  concatenate individual fasta into one big reference sequence file
*/
process combine_refseq_fastas {
  publishDir "${params.refseq_fasta_dir}", mode: 'link'

  label 'process_low'

  input:
  path (individual_fastas) 

  output:
  path ("reference_sequences.fasta"), emit: fasta

  // makes fasta formatted records for each targets
  script:
  """
  cat $individual_fastas > reference_sequences.fasta
  """
}

