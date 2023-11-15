workflow BLAST_UNASSIGNED_SEQUENCES {

  take: 
  unassigned_sequences
  blast_db_dir
  blast_tax_dir
  R_lib_dir

  main:

  ch_versions = Channel.empty()                           

  BLASTN_UNASSIGNED_SEQUENCES(unassigned_sequences, blast_db_dir, blast_tax_dir) 
  ch_versions = ch_versions.mix(BLASTN_UNASSIGNED_SEQUENCES.out.versions)

  ASSIGN_UNASSIGNED_SEQUENCES(BLASTN_UNASSIGNED_SEQUENCES.out.blast_out, unassigned_sequences, R_lib_dir)
  ch_versions = ch_versions.mix(ASSIGN_UNASSIGNED_SEQUENCES.out.versions)

  emit:
  report     = ASSIGN_UNASSIGNED_SEQUENCES.out.unassigned_sequences_report
  versions   = ch_versions
}


/* 
 This process BLASTs all of the unassigned unique sequences 
 (amplicon sequence variants) from dada2 against the NCBI
 nt database to try to figure out what they are
*/

process BLASTN_UNASSIGNED_SEQUENCES {
  publishDir "${params.blast_outdir}", mode: 'link'
  tag "all"

  label 'process_medium'
  label 'error_retry'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  when: 
  params.blast_unassigned_sequences

  input:
  path(unassigned_sequences) 
  path(local_nt_database_dir)  , stageAs: 'local_nt_db_dir'  
  path(blast_tax_dir)          , stageAs: 'local_blast_tax_dir' 

  output:
  path("${unassigned_sequences}.bn_nt"), emit: blast_out 
  path("${unassigned_sequences}")      , emit: unassigned_sequences
  path "versions.yml"                  , emit: versions                                         

  script:                                                                       
  // if BLASTing locally: the location/name of the blast db
  def local_nt_database = "${local_nt_database_dir}/${params.local_nt_database_name}"

  // columns to include in blast output: standard ones plus taxonomy info
  def blastn_columns = "qaccver saccver pident length mismatch gaps qstart qend sstart send evalue bitscore staxid ssciname scomname sblastname sskingdom"
/*
   	    staxid means Subject Taxonomy ID
   	  ssciname means Subject Scientific Name
   	  scomname means Subject Common Name
   	sblastname means Subject Blast Name
   	 sskingdom means Subject Super Kingdom
*/

  def blast_db_params = ""
  if (params.remote_blast_nt) {
    // remote blastn: slower but doesn't require locally installed nt database
    blast_db_params = "-db nt -remote"
  }
  else {
    // local blastn: faster but requires locally installed nt database
    blast_db_params = "-db $local_nt_database"
  }


  """
  #export BLASTDB="$blast_tax_dir"

  # run blastn
  blastn $blast_db_params -task megablast -perc_identity ${params.blast_perc_identity} -qcov_hsp_perc ${params.blast_qcov_hsp_perc} -evalue ${params.max_blast_nt_evalue} -query $unassigned_sequences -outfmt "6 $blastn_columns" -out ${unassigned_sequences}.bn_nt.no_header

  # prepend blast output with the column names so we don't have to manually name them later
  # the perl inline command here is to replace spaces with tabs
  echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header 

  # merge custom header line with actual blast output
  cat blast_header ${unassigned_sequences}.bn_nt.no_header > ${unassigned_sequences}.bn_nt

  cat <<-END_VERSIONS > versions.yml                                          
  "${task.process}":                                                          
      blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')      
  END_VERSIONS          

  """
}

/*
   This process parses the blast output from unassigned sequences and
   populates a spreadsheet with the info. 
*/

process ASSIGN_UNASSIGNED_SEQUENCES {
  label 'process_low'
  tag "all"

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  input:
  path(blast_output) 
  path(unassigned_sequences) 
  path(R_lib_dir_input)                                                         

  output:
  path("non_reference_sequence_assignments.xlsx")  , emit: unassigned_sequences_report
  path "versions.yml"                              , emit: versions                                         

  script:                                                                       
  // only use R lib dir for singularity
  def r_lib_dir = workflow.containerEngine == 'singularity' ? "${R_lib_dir_input}" : "NA"
  """                                                                           
  Rscript ${params.script_dir}/assign_non_ref_seqs.R $r_lib_dir $unassigned_sequences $blast_output
  """             
}
