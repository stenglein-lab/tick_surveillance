#!/usr/bin/env nextflow

/*
    CDC Tick Surveillance Amplicon Sequencing Analysis Pipeline

    December 13, 2021

    Mark Stenglein
*/


// TODO: command line options and parameter checking

// running the pipeline with these params set will output usage info
// e.g. nextflow run main.nf --h
params.help = false
params.h = false

/*
   Pipeline usage output message.

   This will provide usage information in case of a parameter failure
*/
def usageMessage() {

  log.info """

  # Tick-borne pathogen surveillance bioinformatics pipeline

  For more information on this pipeline, see:

   https://github.com/stenglein-lab/tick_surveillance/tree/master/documentation

  ## Pipeline usage:
  
  The typical command for running the pipeline is as follows:

    nextflow run tick_pipeline.nf 
  
  To resume the pipeline (if it stops for some reason or you modify the pipeline scripts): 

    nextflow run tick_pipeline.nf -resume


  ## Pipeline inputs and assumptions:

  For a description of this pipeline's assumptions, input formats, see:

    1. Paired end fastq (or compressed fastq) files in the input/fastq folder.  
       These files should have filenames that end in .fastq or .fastq.gz.
       These files name's should contain the text R1 and R2, corresponding 
       to the 2 paired reads. 

    2. A sample metadata file 


   To document:

    - Expected type of fastq input 
    - Format of primer input 
    - Primer sequences and their orientation.
    - Format of target reference sequences (species name)...
    - Internal control reference sequences
    - Known off target sequences
    - Known off target sequences
    - Versioning of references?
    - NT database info (use RefSeq instead of nt?)
  
  Pipeline parameters: 

  All of the following parameters have default values that can be optionally
  overridden at run time by including the parameter on the command line. For
  example:

   nextflow run tick_pipeline.nf --input_dir a_different_input_directory


    Pipeline input:

    --input_dir                    Input directory. 
                                   [default: ${params.input_dir}] 

    --fastq_dir                    Input directory for fastq files. 
                                   [default: ${params.fastq_dir}] 

    --script_dir                   Directory containing auxiliary pipeline scripts.
                                   [default: ${params.script_dir}] 

    --refseq_dir                   Directory containing target reference sequences.
                                   [default: ${params.refseq_dir}] 

    --targets                      File containing information about the target
                                   and internal control sequences in tsv 
                                   (tab-delimited) format.
                                   [default: ${params.targets}] 

    --primers                      File containing primers used to amplify 
                                   surveillance targets in tsv format.
                                   [default: ${params.primers}] 

    --metadata                     File containing sample metadata
                                   in Excel format.
                                   [default: ${params.metadata}] 


    Pipeline output:

    --outdir                       Output directory into which to place 
                                   result files 
                                   [default: ${params.outdir}] 

    --initial_fastqc_dir           Output directory into which to place 
                                   initial (pre quality-trimming) fastqc 
                                   report files
                                   [default: ${params.initial_fastqc_dir}] 

    --post_trim_fastqc_dir         Output directory into which to place 
                                   post-trimming fastqc report files
                                   [default: ${params.post_trim_fastqc_dir}] 

    --trimmed_outdir               Output directory into which to place 
                                   post-trimming fastqc files
                                   [default: ${params.trimmed_outdir}] 


    Configurable parameters for trimming and assignment:

    --post_trim_min_length         Reads shorter than this after primer trimming
                                   will be discarded from further analysis.
                                   [default: ${params.post_trim_min_length}] 

    --max_blast_refseq_evalue      The maximum blastn E-value for observed
                                   sequences to be considered for possible
                                   assignment to a reference sequence.  
                                   Having an blastn E-value below this is 
                                   necessary but not sufficient to be assigned.
                                   [default: ${params.max_blast_refseq_evalue}] 

    --max_blast_non_refseq_evalue  The maximum blastn E-value for non-reference
                                   sequences to be considered for possible
                                   assignment to a Genbank sequence.  
                                   Having an blastn E-value below this is 
                                   necessary but not sufficient to be assigned.
                                   [default: ${params.max_blast_non_refseq_evalue}] 

    Pipeline help and usage information:

    --help                         Output this usage statement and terminate.
    --h                            Output this usage statement and terminate. 

        """

  // do this extra call to log.info because groovy strips off trailing newlines 
  // from usage message, which makes for slightly less nice looking output.
  log.info """
  """
}

/*

   unimplemented options

    --off_target_refseq_fasta      [NOT IMPLEMENTED]
                                   File containing known off-target amplicon
                                   sequences in fasta format.  Observed sequences
                                   matching these will be removed from analysis.
                                   [default: ${params.off_target_refseq_fasta}] 
*/

/* 
  Pipeline informational output message
*/
// This 
def infoMessage() {

  // tie this versioning to github tags/releases

  log.info """

    Pipeline version:             ${params.pipeline_version}
    Reference sequences version:  ${params.refseq_version}
    Primers version:              ${params.primers_version}

  """
}

/*
  Print usage information if necessary and exit.
*/
if (params.help || params.h) {
    usageMessage()
    exit 0
}


/* 
  Check input parameters 
*/
def check_params_and_input () {

  // check_input_fastq()
  // check_reference_sequences()
  // check_primers()
  // check_metadata()

/*
  if (!$params.input_dir.exists()) {
    log.info """
      Error: input directory $params.input_dir does not exist
    """
    usageMessage()
  }
*/

}


/* 
  Check parameters and input
*/
check_params_and_input()



/* 
  Check for existence of metadata file
*/
Channel
    .fromPath("${params.metadata}", 
                   checkIfExists: true) 
    .set {post_metadata_check_ch}

                                                                                
/* 
  Read in the targets tsv-format file that describes the expected target sequences. 
*/
Channel
    .fromPath(params.targets, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { targets_ch }
                                                                                


/* 
  generate one fasta for each reference sequence
*/
process generate_refseq_fastas {                                                      

  input:                                                                        
  val targets from targets_ch

  output:                                                                        
  path ("${targets.ref_sequence_name}.fasta") into refseq_fastas_ch
  tuple path ("${targets.ref_sequence_name}.fasta"), val (targets) into refseq_fastas_simulate_ch

  // makes fasta formatted records for each target
  script:                                                                       
  """
  printf ">%s\n%s\n" $targets.ref_sequence_name $targets.sequence > ${targets.ref_sequence_name}.fasta
  """
}

/* 
  concatenate individual fasta into one big reference sequence file 
*/
process combine_refseq_fasta {                                                      
  publishDir "${params.refseq_dir}", mode: 'link'                                   

  input:                                                                        
  path (individual_fastas) from refseq_fastas_ch.collect()

  output:                                                                        
  path ("reference_sequences.fasta") into refseq_fasta_ch

  // makes fasta formatted records for each targets 
  script:                                                                       
  """
  cat $individual_fastas > reference_sequences.fasta
  """
}

/*
   Setup indexes and dictionaries needed by downstream processes.

   Only do this once at beginning.
*/
process setup_indexes {
  publishDir "${params.refseq_dir}", mode: 'link'                                   

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::blast=2.12.*" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  input:
  path (refseq_fasta) from refseq_fasta_ch

  output:
  // this output will be a signal that indexes are setup and processes that
  // need them can proceed
  val ("indexes_complete") into post_index_setup_ch
  
  // blast database
  // have to pass paths to both the blast db prefix/name and paths to the additional 
  // files composing the blast db so they will be available to downstream process
  tuple path ("$refseq_fasta"), path ("${refseq_fasta}.*") into refseq_blast_db_ch


  script:
  """
  # Blast database of the reference sequences
  makeblastdb -dbtype nucl -in ${refseq_fasta} -out "${refseq_fasta}"
  
  # make a directory to put the simulated datasets
  mkdir -p ${params.simulated_fastq_dir}
  """
}


// the sizes of control datasets to make
// TODO: how to make this parameterized but not with a list type 
// control_dataset_sizes_ch  = Channel.of(params.simulated_dataset_sizes) 
control_dataset_sizes_ch  = Channel.of(0, 100, 1000)

/* 
  generate simulated fastq reads for each reference sequence 
  this will serve as an internal control
*/
process simulate_refseq_fastq {
  publishDir "${params.simulated_fastq_dir}", mode: 'link'                                   

  input:                                                                        
  tuple val (dataset_size), path(ref_fasta), val(target) from control_dataset_sizes_ch.combine(refseq_fastas_simulate_ch)

  output:                                                                        
  tuple val ("${target.ref_sequence_name}_${dataset_size}"), path ("${target.ref_sequence_name}_${dataset_size}*.fastq") into simulated_fastq_ch

  when:
  params.use_simulated_test_datasets

  script:                                                                       

  // this is how the simulated_fastq will be named
  def fastq_prefix =  "${target.ref_sequence_name}_${dataset_size}"
  
  """
  ${params.script_dir}/simulate_fastq.pl -n $dataset_size -f $ref_fasta -pre $fastq_prefix -l ${params.simulated_read_length} -e ${params.simulated_error_profile_file}
  # ${params.script_dir}/simulate_fastq.pl -n $dataset_size -f $ref_fasta -pre $fastq_prefix -l ${params.simulated_read_length} 
  """
}

// TODO: signal to proceed once datasets generated

/*
 These fastq files represent the main input to this workflow
 
 Expecting files with _R1 or _R2 in their names corresponding to paired-end reads
*/

Channel
    .fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq*", 
                   size: 2, 
                   // checkIfExists: true, 
                   maxDepth: 1)
    .into {samples_ch_qc; samples_ch_trim}


/* 

 This channel generates the primer sequences that were 
 used to amplify surveillance targets.  These primer sequences will be trimmed 
 off of read pairs and used to identify legitimate PCR products, which will
 contain an expected F/R primer pair at the ends. 

*/
Channel
    .fromPath(params.primers, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { primers_ch }


/*
   this combinatorially mixes the fastq file pairs with the primer sequences
   to be trimmed, creating a new channel with all possible combinations of
   input fastq and primer pair
*/

// primers_and_samples = primers_ch.combine(samples_ch_trim)
// primers_and_samples = primers_ch.combine(samples_ch_trim.concat(simulated_fastq_ch))

primers_and_samples = params.use_simulated_test_datasets ?  primers_ch.combine(simulated_fastq_ch) : primers_ch.combine(samples_ch_trim)

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

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
  } else {
      container "quay.io/biocontainers/fastqc:0.11.9--0"
  }

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  val(sample_id) into post_initial_qc_ch

  script:
  """
  which fastqc
  mkdir -p  ${params.initial_fastqc_dir} 
  fastqc -o ${params.initial_fastqc_dir} $initial_fastq 
  """
}

/*
 Use multiqc to merge initial fastqc reports

 TODO: will this report old, accumulated fastqc reports if the pipeline is re-run without cleaning up the work directory?
*/
process initial_multiqc {
  publishDir "${params.outdir}", mode:'link'

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::multiqc=1.11" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
  }

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

process trim_primer_seqs {                                                      
  label 'lowmem'

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::cutadapt=3.5" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
  }
                                                                              
  input:                                                                      
  val("indexes_complete") from post_index_setup_ch
  tuple val(primers), val(sample_id), path(initial_fastq) from primers_and_samples

  output:                                                                     
  tuple val(sample_id), path("*.R1_${primers.primer_name}.fastq.gz"), path("*.R2_${primers.primer_name}.fastq.gz") into primer_trimmed_ch_ungrouped
                                                                              
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
  cutadapt $primer_args --discard-untrimmed  -e ${params.amplicon_primers_max_error_fraction} --minimum-length ${params.post_trim_min_length} $f1 $f2 -o ${sample_id}.R1_${primers.primer_name}.fastq.gz -p ${sample_id}.R2_${primers.primer_name}.fastq.gz
  """                                                                         
}
                                                                                
/*
  This process concatenates the outputs from the individual cutadapt calls

  Cutadapt is run multiple times per fastq pair to create individual outputs 
  containing read pairs that contained an expected primer pair at their ends

  This step concatenates all those outputs into a single file per original
  fastq pair. 

  Because *all* cutadapt output flows into the primer_trimmed_ch_ungrouped channel,
  it is necessary to bin these outputs into groups corresponding to the original
  fastq read pairs.  

  The trick to doing this is the .groupTuple() operator, which
  groups tuples (sample_id + read1_trimmed_fastq + read2_trimmed_fastq) based on 
  a shared sample_id.

  This process also does one more round of trimming to remove stanard Illumina
  adapter sequences

*/
process collect_cutadapt_output {                                               
  publishDir "${params.trimmed_outdir}", mode:'link'                                    
                                                                                
  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::cutadapt=3.5" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
  }

  input:
  // the groupTuple() operator here will consolidate tuples with a shared sample_id
  tuple val(sample_id), path(individual_r1), path(individual_r2) from primer_trimmed_ch_ungrouped.groupTuple()
                                                                                
  output:                                                                       
  val(sample_id) into post_trim_ch
  val(sample_id) into post_trim_sample_ch
  tuple val(sample_id), path("*_trimmed.fastq.gz") into post_trim_qc_ch
                                                                                
  script:                                                                       

  // cutadapt parameters to trim TruSeq-style adapters
  // see: https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage
  def truseq_cutadapt = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  
  // this shortened version of the Nextera adapter found in some amplicons
  def nextera_cutadapt = "-a CTGTCTCTTATACA -A CTGTCTCTTATACA"

  // these will be the names of the new concatenated outputs, before truseq trimming
  def f1 = "${sample_id}.R1_individual_primers.fastq.gz"
  def f2 = "${sample_id}.R2_individual_primers.fastq.gz"

  """                                                                           
  # concatenate the individually-trimmed files (one for each amplicon target)
  # note that cat will concatenate .gz files just fine
  cat $individual_r1 > $f1
  cat $individual_r2 > $f2

  # one last cutadapt command to remove adapter sequences and too-short read pairs
  cutadapt $nextera_cutadapt \
    -O ${params.adapters_min_overlap} \
    --minimum-length ${params.post_trim_min_length} \
    -o ${sample_id}_R1_trimmed.fastq.gz \
    -p ${sample_id}_R2_trimmed.fastq.gz \
    ${f1} ${f2}
  """                                                                           
}         


/*
 Use fastqc to do QC on post-trimmed fastq
*/
process post_trim_qc {
  label 'lowmem'

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
  } else {
      container "quay.io/biocontainers/fastqc:0.11.9--0"
  }

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

 This report should show, for instance, the absence of any remaining
 Illumina adapter sequences. 

 TODO: will this report old fastqc reports if the pipeline is re-run?

*/
process post_trim_multiqc {
  publishDir "${params.outdir}", mode: 'link'

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::multiqc=1.11" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
  }

  input:
  val(all_sample_ids) from post_trim_multiqc_ch.collect()

  output:
  path("post_trim_qc_report.html")

  """
  multiqc -n "post_trim_qc_report.html" -m fastqc ${params.post_trim_fastqc_dir}
  """
}


/* 
  Run dada2 on trimmed read pairs to:

 - perform error-correction on misscalled bases
 - merge read pairs
 - cluster identical read pairs
 - remove chimeric sequences
 - tabulate the frequencies of particular unique sequences

*/

process run_dada_on_trimmed {
  publishDir "${params.outdir}", mode: 'link'

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22.*" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bioconductor-dada2%3A1.22.0--r41h399db7b_0"
  } else {
      container "quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0"
  }

  input:
  val(all_sample_ids) from post_trim_ch.collect()

  output:
  path("dada_seqtab.txt") into post_dada_run_ch

  script:                                                                       
  """                                                                             
  # Run dada2 using trimmed fastq as input and create a tabular output of results
  Rscript ${params.script_dir}/run_dada_on_trimmed.R ${params.script_dir} ${params.trimmed_outdir}
  """             
}


/* 
  Tidy up dada2 output and write out observed sequences (ASVs) in fasta format
*/

process tidy_dada_output {
  publishDir "${params.outdir}", mode: 'link'

  // conda / singularity info for this process
  conda (params.enable_conda ? "conda-forge::r-base=4.1.* conda-forge::r-tidyverse=1.3.* conda-forge::r-openxlsx=4.2.*" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "library://stenglein-lab/r_tools/r_tools:1.0.0"
  } else {                                                                      
      container "library://stenglein-lab/r_tools/r_tools:1.0.0"
  }        

  input:
  path(dada_output) from post_dada_run_ch

  output:
  tuple path("observed_sequences.fasta"), path("sequence_abundance_table.tsv") into post_dada_ch

  script:                                                                       
  """                                                                             
  # This R script creates an output named observed_sequences.fasta containing all of the
  # unique sequences observed in the amplicon dataset
  # and a sequence_abundance_table.tsv, which lists the abundances of these
  # sequences in each dataset
  Rscript ${params.script_dir}/tidy_dada_output.R ${params.script_dir} $dada_output
  """             
}

/* 
 This process aligns all of the unique sequences reported by dada2
 to the set of expected reference sequences using blastn.
*/
process compare_observed_sequences_to_ref_seqs {
  publishDir "${params.outdir}", mode: 'link'

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::blast=2.12.*" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  input:
  tuple path(sequences), path(abundance_table) from post_dada_ch
  tuple path(refseq_blast_db), path (refseq_blast_db_files) from refseq_blast_db_ch

  output:
  tuple path(abundance_table), path("${sequences}.bn_refseq") into post_compare_ch

  script:                                                                       
  // this is similar to the default blastn output except gaps replaces gapopens, because seems more useful!
  // also include other useful columns like query length, subject length, etc
  def blastn_columns = "qaccver saccver pident length mismatch gaps qlen slen bitscore"
  """                                                                           
  blastn -db ${refseq_blast_db} -task blastn -evalue ${params.max_blast_refseq_evalue} -query $sequences -outfmt "6 $blastn_columns" -out ${sequences}.bn_refseq.no_header
  # prepend blast output with the column names so we don't have to manually name them later
  echo $blastn_columns > blast_header.no_perl
  echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header 
  cat blast_header ${sequences}.bn_refseq.no_header > ${sequences}.bn_refseq
  """             
}

/*
  This process takes the blast alignment information from the above process
  and decides whether those sequence are sufficiently similar to the reference 
  sequences to be assigned to them.
*/
process assign_observed_sequences_to_ref_seqs {
  publishDir "${params.outdir}", mode: 'link'


  // conda / singularity info for this process
  conda (params.enable_conda ? "conda-forge::r-base=4.1.* conda-forge::r-tidyverse=1.3.* conda-forge::r-openxlsx=4.2.*" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "library://stenglein-lab/r_tools/r_tools:1.0.0"
  } else {                                                                      
      container "library://stenglein-lab/r_tools/r_tools:1.0.0"
  }        
  

  input:
  path(metadata) from post_metadata_check_ch
  tuple path(abundance_table), path(blast_output) from post_compare_ch

  output:
  path("unassigned_sequences.fasta") into post_assign_to_refseq_ch
  path("identified_targets.xlsx")
  path("identified_targets.tsv")

  script:                                                                       
  """                                                                           
  Rscript ${params.script_dir}/assign_observed_seqs_to_ref_seqs.R ${params.script_dir} $abundance_table $blast_output $metadata ${params.targets} 
  """             
}


/* 
 This process BLASTs all of the unassigned unique sequences 
 (amplicon sequence variants) from dada2 against the NCBI
 nt database to try to figure out what they are
*/
process blast_unassigned_sequences {
  publishDir "${params.outdir}", mode: 'link'

  // conda / singularity info for this process
  conda (params.enable_conda ? "bioconda::blast=2.12.*" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  input:
  path(unassigned_sequences) from post_assign_to_refseq_ch

  output:
  path("${unassigned_sequences}.bn_nt") into post_blast_unassigned_ch
  path("${unassigned_sequences}") into post_blast_unassigned_seq_ch

  script:                                                                       
  // columns to include in blast output: standard ones plus taxonomy info
  def blastn_columns = "qaccver saccver pident length mismatch gaps qstart qend sstart send evalue bitscore staxid ssciname scomname sblastname sskingdom"
/*
   	    staxid means Subject Taxonomy ID
   	  ssciname means Subject Scientific Name
   	  scomname means Subject Common Name
   	sblastname means Subject Blast Name
   	 sskingdom means Subject Super Kingdom
*/

  if (params.remote_blast_nt) {

    // remote blastn: slower but doesn't require locally installed nt database

    """

    # run remote blastn
    blastn -db nt -task megablast -remote -evalue ${params.max_blast_nt_evalue} -query $unassigned_sequences -outfmt "6 $blastn_columns" -out ${unassigned_sequences}.bn_nt.no_header

    # prepend blast output with the column names so we don't have to manually name them later
    # the perl inline command here is to replace spaces with tabs
    echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header 

    # merge custom header line with actual blast output
    cat blast_header ${unassigned_sequences}.bn_nt.no_header > ${unassigned_sequences}.bn_nt

    """

  }
  else {

    // local blastn: faster but requires locally installed nt database

    """

    # run local blastn
    blastn -db nt -task megablast -evalue ${params.max_blast_nt_evalue} -query $unassigned_sequences -outfmt "6 $blastn_columns" -out ${unassigned_sequences}.bn_nt.no_header

    # prepend blast output with the column names so we don't have to manually name them later
    # the perl inline command here is to replace spaces with tabs
    echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header 

    # merge custom header line with actual blast output
    cat blast_header ${unassigned_sequences}.bn_nt.no_header > ${unassigned_sequences}.bn_nt

    """             
  }
}

/*
   This process parses the blast output from unassigned sequences and
   populates a spreadsheet with the info. 
*/
process assign_non_ref_seqs {
  publishDir "${params.outdir}", mode: 'link'

  // conda / singularity info for this process
  conda (params.enable_conda ? "conda-forge::r-base=4.1.* conda-forge::r-tidyverse=1.3.* conda-foge::r-openxlsx=4.2.*" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "library://stenglein-lab/r_tools/r_tools:1.0.0"
  } else {                                                                      
      container "library://stenglein-lab/r_tools/r_tools:1.0.0"
  }        

  input:
  path(blast_output) from post_blast_unassigned_ch
  path(unassigned_sequences) from post_blast_unassigned_seq_ch

  output:
  path("non_reference_sequence_assignments.xlsx")

  script:                                                                       
  """                                                                           
  Rscript ${params.script_dir}/assign_non_ref_seqs.R ${params.script_dir} $unassigned_sequences $blast_output
  """             
}

