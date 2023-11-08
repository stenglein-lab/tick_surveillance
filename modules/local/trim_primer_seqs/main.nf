/*
   A function to return the complement of a DNA sequence
   combine with reverse() to get the reverse complement

   from http://groovyconsole.appspot.com/script/29005
*/
String.metaClass.complement = {
  def complements = [ A:'T', T:'A', U:'A', G:'C',
                      C:'G', Y:'R', R:'Y', S:'S',
                      W:'W', K:'M', M:'K', B:'V',
                      D:'H', H:'D', V:'B', N:'N' ]
  delegate.toUpperCase().collect { base ->
    complements[ base ] ?: 'X'
  }.join()
}

/*
 Trim primer sequences.

 This process will look output read pairs that have one of the expected primer pairs at
 the expected positions at the ends of R1 and R2.  These primers will be trimmed and
 matching read pairs will be output to a new file.  This avoid amplification products that
 were formed by unexpected primer combinations (since this is a multiplex PCR assay).

 This trimming uses cutadapt as described in:

 https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads

 This step also throws away read pairs shorter than ${params.post_trim_min_length}

 This process will be run once per input fastq file pair per primer pair

*/

process TRIM_PRIMER_SEQS {
  label 'many_forks'
  label 'process_low'
  tag   "${meta.id}/${primers.primer_name}"
  publishDir "${params.cutadapt_trim_reports}", pattern: '*_summary.txt', mode: 'copy'

  // if using conda
  conda "bioconda::cutadapt=3.5"

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
  }
                       
  input:
  tuple val(primers), val(meta), path(initial_fastq) 
 
  output:                                                                     
  // output the R1 and R2s as seperate list elements rather than a combined sub-list
  // because we want to handle them separately later
  tuple val(meta), 
        path("*R1_${primers.primer_name}.fastq.gz"), 
        path("*R2_${primers.primer_name}.fastq.gz")         , emit: trimmed_reads
  path "versions.yml"                                       , emit: versions
  path("${meta.id}_${primers.primer_name}_summary.txt")     , emit: cutadapt_trim_report                                         

                                             
  script:
  def primer_f_rc = primers.primer_f_seq.reverse().complement()
  def primer_r_rc = primers.primer_r_seq.reverse().complement()
  def f1 = initial_fastq[0]
  def f2 = initial_fastq[1]

  // setup arguments for trimming this particular pair of primers
  // see: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads
  def primer_args = "-a " + primers.primer_f_seq + "..." + primer_r_rc + " -A " + primers.primer_r_seq + "..."  + primer_f_rc


  """
  # the --discard-untrimmed option means that the output of this command will contain
  # only those reads that match this particular expected primer pair because only
  # those matching read pairs will be trimmed, and non-trimmed pairs will get discarded.
  #
  # from the cutadapt manual:
  # "Use --discard-untrimmed to throw away all read pairs in which R1 
  # doesn’t start with FWDPRIMER or in which R2 does not start with REVPRIMER"
  #
  cutadapt $primer_args \
    --discard-untrimmed  \
    -e ${params.amplicon_primers_max_error_fraction} \
    --minimum-length ${params.post_trim_min_length} \
    $f1 \
    $f2 \
    -o ${meta.id}.R1_${primers.primer_name}.fastq.gz \
    -p ${meta.id}.R2_${primers.primer_name}.fastq.gz \
    > ${meta.id}_${primers.primer_name}_summary.txt

  cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
    cutadapt: \$(cutadapt --version)
  END_VERSIONS         
  """
}
