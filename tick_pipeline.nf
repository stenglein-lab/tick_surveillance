#!/usr/bin/env nextflow

/*
    CDC Tick Surveillance Amplicon Sequencing Analysis Pipeline

    October 22, 2020 

    Mark Stenglein
*/


// TODO: command line options and parameter checking


// --------------------------------
// Directories for input and output
// --------------------------------
params.fastq_dir = "$baseDir/a_few_fastq/"
params.initial_fastqc_dir = "$baseDir/initial_fastqc/" 
params.post_trim_fastqc_dir = "$baseDir/post_trim_fastqc/" 
params.outdir = "$baseDir/results"                                                       
params.trimmed_outdir = "$baseDir/trimmed_fastq"                                                       

// where are R scripts found...                                                 
params.R_bindir="${baseDir}/scripts"  

// -------------------
// Reference sequences
// -------------------
params.refseq_dir = "refseq"
params.refseq_fasta = "${params.refseq_dir}/reference_sequences.fasta"
params.refseq_version = "2020-11-11"

params.primers = "${params.refseq_dir}/primers.csv"


// ------------------
// Trimming settings
// ------------------
// shortest amplicon = tick Actin @ 196 bp
params.post_trim_min_length = "100" 




// TODO: command line arg processing and validating 

// TODO: check that appropriate refseq files exist (fasta, gb, etc.)

// TODO: handle single end or paired end, possibly automatically

// TODO: better control of concurrency in terms of processes

// TODO: move scripts to a bin subdir

// TODO: put all up on github

// TODO: conda

// TODO: BSQR fails in case of no mapping reads... deal with this possibility 

/*
 These fastq files represent the main input to this workflow
 
 Expecting files with _R1 or _R2 in their names corresponding to paired-end reads
*/
Channel
    .fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq", size: 2, checkIfExists: true, maxDepth: 1)
    .into {samples_ch_qc; samples_ch_trim}

/* 
 The primer sequences that will be trimmed off of read pairs

 Only read pairs that have an expected matching pair of primer sequences
 will be retained.
*/
Channel                                                                         
    .fromPath(params.primers)                                                   
    .splitCsv(header:true)                                                      
    .map{ row-> tuple(row.target, row.primer_r_name, row.primer_f_seq, row.primer_r_name, row.primer_r_seq) }
    .set { primers_ch }   
                                                                                
// this combinatorially combines the fastq file pairs with the primer sequences
// to be trimmed, creating a new channel
primers_and_samples = primers_ch.combine(samples_ch_trim)   

/*
   Setup some initial indexes and dictionaries needed by downstream processes.
   Only do this once at beginning.
   TODO: remove initial .dict file
*/
/*
process setup_indexes {

  output:
  // this output will be a signal that indexes are setup and processes that
  // need them can proceed
  val("indexes_complete") into post_index_setup_ch

  script:
  """
  # -----------------------
  # bwa index viral refseq
  # -----------------------
  bwa index ${params.viral_fasta}

  # -----------------
  # GATK index setup
  # -----------------
  # setup gatk indexes for BSQR
  ${params.gatk_exe} IndexFeatureFile -I ${params.ignore_regions} 

  rm -f "${baseDir}/viral_refseq/${params.viral_refseq_name}.dict"
  ${params.gatk_exe} CreateSequenceDictionary -R ${params.viral_fasta}
  """
}
*/


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
 Run fastqc on input fastq 
*/
process initial_qc {
  label 'lowmem'

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  val(sample_id) into post_initial_qc_ch
  val(sample_id) into write_datasets_ch
  // TODO: count

  script:
  """
  which fastqc
  mkdir -p  ${params.initial_fastqc_dir} 
  fastqc -o ${params.initial_fastqc_dir} $initial_fastq 
  """
}

/*
 Use multiqc to merge initial fastqc reports
*/
process initial_multiqc {
  publishDir "${params.outdir}", mode:'link'

  input:
  val(all_sample_ids) from post_initial_qc_ch.collect()

  output: 
  path("initial_qc_report.html")

  script:
  """
  multiqc -n "initial_qc_report.html" -m fastqc ${params.initial_fastqc_dir}
  """
} 

/*
 trim primer sequences 
*/
process trim_primer_seqs {                                                      
  label 'lowmem'
                                                                              
  input:                                                                      
  tuple val(target), 
        val(primer_f_name), val(primer_f), 
        val(primer_r_name), val(primer_r), 
        val(sample_id), path(initial_fastq) from primers_and_samples
                                                                              
  output:                                                                     
  tuple val(sample_id), path("*.R1_${target}.fastq"), path("*.R2_${target}.fastq") into primer_trimmed_ch_ungrouped
                                                                              
  script:                                                                     
  def primer_f_rc = primer_f.reverse().complement()                           
  def primer_r_rc = primer_r.reverse().complement()                           
  def f1 = initial_fastq[0]
  def f2 = initial_fastq[1]

  // the cutadapt primer trimming arguments should be of the form: 
  // -a ^FWDPRIMER...RCREVPRIMER -A ^REVPRIMER...RCFWDPRIMER
  // see: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads
  // def primer_args = "-a " + target + "_F=^" + primer_f + "..." + primer_r_rc + " -A " + target + "_R=^" + primer_r + "..."  + primer_f_rc
  def primer_args = "-a " + target + "_F=^" + primer_f_rc + "..." + primer_r + " -A " + target + "_R=^" + primer_r_rc + "..."  + primer_f
                                                                              
  """                                                                         
  # the --discard-untrimmed option means that the output of this command will contain
  # only those reads that match this particular expected primer pair
  #
  # from the cutadapt manual:
  # "Use --discard-untrimmed to throw away all read pairs in which R1 
  # doesn’t start with FWDPRIMER or in which R2 does not start with REVPRIMER"
  #
  cutadapt $primer_args --discard-untrimmed  --minimum-length ${params.post_trim_min_length} $f1 $f2 -o ${sample_id}.R1_${target}.fastq -p ${sample_id}.R2_${target}.fastq
  """                                                                         
}
                                                                                
/*
  This process concatenates the outputs from the individual cutadapt calls
*/
process collect_cutadapt_output {                                               
  publishDir "${params.trimmed_outdir}", mode:'link'                                    
                                                                                
  input:
  // the groupTuple() operator here will consolidate tuples with a shared sample_id
  tuple val(sample_id), path(individual_r1), path(individual_r2) from primer_trimmed_ch_ungrouped.groupTuple()
                                                                                
  output:                                                                       
  // tuple val(sample_id), path("*_trimmed.fastq") into post_trim_ch
  val(sample_id) into post_trim_ch
  tuple val(sample_id), path("*_trimmed.fastq") into post_trim_qc_ch
                                                                                
  script:                                                                       

  // cutadapt parameters to trim TruSeq-style adapters
  // see: https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage
  def truseq_cutadapt = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  def f1 = "${sample_id}.r1_individual_primers.fastq"
  def f2 = "${sample_id}.r2_individual_primers.fastq"

  """                                                                           
  # concatenate the individually-trimmed files (one for each amplicon target)
  cat $individual_r1 > $f1
  cat $individual_r2 > $f2

  # one last cutadapt command to remove adapter sequences and too-short read pairs
  cutadapt $truseq_cutadapt --discard-trimmed --minimum-length ${params.post_trim_min_length} \
    -o ${sample_id}_R1_trimmed.fastq -p ${sample_id}_R2_trimmed.fastq ${f1} ${f2}

  """                                                                           
}         

/*
 Use fastqc to do QC on post-trimmed fastq
*/
process post_trim_qc {
  label 'lowmem'

  input:
  tuple val(sample_id), path(input_fastq) from post_trim_qc_ch

  output:
  val(sample_id) into post_trim_multiqc_ch

  script:

  """
  mkdir -p  ${params.post_trim_fastqc_dir} 
  fastqc -o ${params.post_trim_fastqc_dir} $input_fastq
  """
}

/*
 Use multiqc to merge post-trimming fastq reports
*/
process post_trim_multiqc {
  publishDir "${params.outdir}", mode: 'link'

  input:
  val(all_sample_ids) from post_trim_multiqc_ch.collect()

  output:
  path("post_trim_qc_report.html")

  """
  multiqc -n "post_trim_qc_report.html" -m fastqc ${params.post_trim_fastqc_dir}
  """
}


process run_dada_on_trimmed {
  publishDir "${params.outdir}", mode: 'link'

  input:
  val(all_sample_ids) from post_trim_ch.collect()

  output:
  path("wide_seqtab.txt")

  script:                                                                       
  """                                                                           
  Rscript ${params.R_bindir}/run_dada_on_trimmed.R ${params.R_bindir} ${params.trimmed_outdir}
  """             


}

