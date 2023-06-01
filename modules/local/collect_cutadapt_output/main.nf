/*
  This process concatenates the outputs from the individual cutadapt calls

  Cutadapt is run multiple times per fastq pair to create individual outputs 
  containing read pairs that contained an expected primer pair at their ends

  This step concatenates all those outputs into a single file per original
  fastq pair. 

  This process also does one more round of trimming to remove standard Illumina
  adapter sequences
*/
process COLLECT_CUTADAPT_OUTPUT {      
  label      'process_low'
  tag        "${meta.id}"
  publishDir "${params.trimmed_outdir}", pattern: "*.fastq.gz"
                                                                   
  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
  }

  input:
  tuple val(meta), path(individual_R1_fastq), path(individual_R2_fastq)
                                                                                
  output:                                                                       
  tuple val(meta), path("*_trimmed.fastq.gz"), emit: trimmed_reads 
  path "versions.yml"                        , emit: versions                                         
                                                          
  script:                                                                       

  // cutadapt parameters to trim TruSeq-style adapters
  // see: https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage
  def truseq_cutadapt = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  
  // this shortened version of the Nextera adapter found in some amplicons
  def nextera_cutadapt = "-a CTGTCTCTTATACA -A CTGTCTCTTATACA"

  // these will be the names of the new concatenated outputs, before truseq trimming
  def f1 = "${meta.id}.R1_individual_primers.fastq.gz"
  def f2 = "${meta.id}.R2_individual_primers.fastq.gz"

  """
  # concatenate the individually-trimmed files (one for each amplicon target)
  # note that cat will concatenate .gz files just fine
  cat $individual_R1_fastq > $f1
  cat $individual_R2_fastq > $f2

  # one last cutadapt command to remove adapter sequences and too-short read pairs
  cutadapt $nextera_cutadapt \
    -O ${params.adapters_min_overlap} \
    --minimum-length ${params.post_trim_min_length} \
    -o ${meta.id}_R1_trimmed.fastq.gz \
    -p ${meta.id}_R2_trimmed.fastq.gz \
    ${f1} ${f2}

  cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
    cutadapt: \$(cutadapt --version)
  END_VERSIONS         
  """
}
