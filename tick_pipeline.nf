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
params.input_dir = "$baseDir/input/"
params.fastq_dir = "${params.input_dir}/fastq/"

// Directory where R (and any other) scripts are.                                                 
params.script_dir="${baseDir}/scripts"  

params.outdir = "$baseDir/results"                                                       
params.initial_fastqc_dir = "${params.outdir}/initial_fastqc/" 
params.post_trim_fastqc_dir = "${params.outdir}/post_trim_fastqc/" 
params.trimmed_outdir = "${params.outdir}/trimmed_fastq"                                                       


// the Local installation of NCBI nt database
params.local_nt_database ="/scicomp/reference/ncbi-blast/current/nt"


// -------------------
// Reference sequences
// -------------------
// TODO: check that appropriate refseq files exist
params.refseq_dir = "${baseDir}/refseq"
params.targets = "${params.refseq_dir}/targets.tsv"
params.off_target_refseq_fasta = "${params.refseq_dir}/off_target_products.fasta"

// ---------------
// Sample metadata
// ---------------
// TODO (possibly) Force the user to specify a file (?)
params.metadata = "${params.input_dir}/sample_metadata.xlsx"
// params.metadata = "${params.input_dir}/Group2_metadata.xlsx"

// ------------------------------------
// Simulated internal-control datasets
// ------------------------------------
params.simulated_fastq_dir = "${params.refseq_dir}/simulated_fastq/"

// simulated internal control datasets
// background probability of errors for simulating fastq
params.simulated_error_profile_file = "${baseDir}/refseq/q35_errors.txt"
params.simulated_read_length = 250


// --------
// Versions
// --------
// TODO: versioning of pipeline, reference sequences, etc: how to best do (?)
params.pipeline_version = "2021-05-19"
params.refseq_version = "2021-01-08"
params.primers_version = "2021-01-08"

// ---------
// Trimming 
// ---------
// the primers that will be trimmed off of the ends of amplicons
params.primers = "${params.refseq_dir}/primers.tsv"

// shortest amplicon = tick Actin @ 196 bp
params.post_trim_min_length = "100" 

// to trim Tru-seq style adapters need 10 bp of overlap
// this to avoid false-positive trimming
params.adapters_min_overlap = 10


// ----------------------------------------------
// Blast e-value cutoffs and other configuration
// ----------------------------------------------
params.max_blast_refseq_evalue = "1e-10"
params.max_blast_nt_evalue = "1e-10"

// Should BLAST of unassigned sequences be remote (use -remote option)
// This will cause blastn to run much slower but will not require
// the installation of a local database.
//
// This is the option that should default to unless a local database
// is specified and verified to exist.
params.remote_blast_nt = true

// Run testing or not
params.do_testing = true

// these will output usage info
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

    2. An optional sample metadata file 


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

  // how to best handle this versioning?
  // via github?

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

  Use this to generate internal control datasets. 
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

  // makes fasta formatted records for each targets 
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

  input:
  path (refseq_fasta) from refseq_fasta_ch

  output:
  // this output will be a signal that indexes are setup and processes that
  // need them can proceed
  val ("indexes_complete") into post_index_setup_ch
  
  // blast database
  path ("${refseq_fasta}") into refseq_blast_db_ch

  script:
  """
  # Blast database of the reference sequences
  makeblastdb -dbtype nucl -in ${refseq_fasta} -out "${refseq_fasta}"
  
  # make a directory to put the simulated datasets
  mkdir -p ${params.simulated_fastq_dir}
  """
}


// the sizes of control datasets to make
control_dataset_sizes_ch  = Channel.from(1, 10, 100)

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

  script:                                                                       

  // this is how to the simulated_fastq will be named
  def fastq_prefix =  "${target.ref_sequence_name}_${dataset_size}"
  
  """
  # ${params.script_dir}/simulate_fastq.pl -n $dataset_size -f $ref_fasta -pre $fastq_prefix -l ${params.simulated_read_length} -e ${params.simulated_error_profile_file}
  ${params.script_dir}/simulate_fastq.pl -n $dataset_size -f $ref_fasta -pre $fastq_prefix -l ${params.simulated_read_length} 
  """
}


/*
 These fastq files represent the main input to this workflow
 
 Expecting files with _R1 or _R2 in their names corresponding to paired-end reads
*/

Channel
    .fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq*", 
                   size: 2, 
                   checkIfExists: true, 
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
primers_and_samples = primers_ch.combine(samples_ch_trim.concat(simulated_fastq_ch))   



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
  // TODO: count (?)

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

 This process will be run once per input fastq file pair per read pair

 TODO: There is an assumption here that the primer sequences are entered in the correct 
       orientation relative to the Illumina read1 and read2.  
       Unit testing should include a check for that.

*/

process trim_primer_seqs {                                                      
  label 'lowmem'
                                                                              
  input:                                                                      
  val("indexes_complete") from post_index_setup_ch
  tuple val(primers), val(sample_id), path(initial_fastq) from primers_and_samples

/*
  tuple val(primer_name), 
        val(primer_f_name), val(primer_f), 
        val(primer_r_name), val(primer_r), 
        val(sample_id), path(initial_fastq) from primers_and_samples
*/
                                                                              
  output:                                                                     
  tuple val(sample_id), path("*.R1_${primers.primer_name}.fastq.gz"), path("*.R2_${primers.primer_name}.fastq.gz") into primer_trimmed_ch_ungrouped
                                                                              
  script:                                                                     
  def primer_f_rc = primers.primer_f_seq.reverse().complement()                           
  def primer_r_rc = primers.primer_r_seq.reverse().complement()                           
  def f1 = initial_fastq[0]
  def f2 = initial_fastq[1]

  // setup arguments for trimming this particular pair of primers
  // see: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads
  
  // def primer_args = "-a " + primers.primer_name + "_F=^" + primer_f_rc + "..." + primers.primer_r_seq + " -A " + primers.primer_name  + "_R=^" + primer_r_rc + "..."  + primers.primer_f_seq
  // def primer_args = "-a " + primers.primer_name + "_F=" + primer_f_rc + "..." + primers.primer_r_seq + " -A " + primers.primer_name  + "_R=" + primer_r_rc + "..."  + primers.primer_f_seq
  def primer_args = "-a " + primers.primer_name + "_F=" + primers.primer_f_seq + "..." + primer_r_rc + " -A " + primers.primer_name  + "_R=" + primers.primer_r_seq + "..."  + primer_f_rc
                                                                              
  """                                                                         
  # the --discard-untrimmed option means that the output of this command will contain
  # only those reads that match this particular expected primer pair because only
  # those matching read pairs will be trimmed, and non-trimmed pairs will get discarded.
  #
  # from the cutadapt manual:
  # "Use --discard-untrimmed to throw away all read pairs in which R1 
  # doesn’t start with FWDPRIMER or in which R2 does not start with REVPRIMER"
  #
  cutadapt $primer_args --discard-untrimmed  --minimum-length ${params.post_trim_min_length} $f1 $f2 -o ${sample_id}.R1_${primers.primer_name}.fastq.gz -p ${sample_id}.R2_${primers.primer_name}.fastq.gz
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
                                                                                
  input:
  // the groupTuple() operator here will consolidate tuples with a shared sample_id
  tuple val(sample_id), path(individual_r1), path(individual_r2) from primer_trimmed_ch_ungrouped.groupTuple()
                                                                                
  output:                                                                       
  val(sample_id) into post_trim_ch
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

  input:
  val(all_sample_ids) from post_trim_ch.collect()

  output:
  path("observed_sequences.fasta") into post_dada_seq_ch
  path("sequence_abundance_table.tsv") into post_dada_tidy_ch

  script:                                                                       
  """                                                                             
  # This R file creates an output named observed_sequences.fasta containing all of the
  # unique sequences observed in the amplicon dataset
  # and a tidy_sequence_table.tsv, which lists the abundances of these
  # sequences in each dataset
  Rscript ${params.script_dir}/run_dada_on_trimmed.R ${params.script_dir} ${params.trimmed_outdir}
  """             
}

/* 
 This process aligns all of the unique sequences reported by dada2
 to the set of expected reference sequences using blastn.
*/
process compare_observed_sequences_to_ref_seqs {
  publishDir "${params.outdir}", mode: 'link'

  input:
  path(sequences) from post_dada_seq_ch
  path(tidy_table) from post_dada_tidy_ch
  path(refseq_blast_db) from refseq_blast_db_ch

  output:
  path("${sequences}.bn_refseq") into post_compare_ch

  script:                                                                       
  // this is almost the default blastn output except gaps replaces gapopens, because seems more useful!
  def blastn_columns = "qaccver saccver pident length mismatch gaps qstart qend sstart send evalue bitscore"
  """                                                                           
  blastn -db ${params.refseq_dir}/${refseq_blast_db} -task blastn -evalue ${params.max_blast_refseq_evalue} -query $sequences -outfmt "6 $blastn_columns" -out ${sequences}.bn_refseq.no_header
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

  input:
  path(metadata) from post_metadata_check_ch
  path(blast_output) from post_compare_ch
  path(tidy_table) from post_dada_tidy_ch

  output:
  path("unassigned_sequences.fasta") into post_assign_to_refseq_ch
  path("identified_targets.xlsx")
  path("identified_targets.tsv")

  script:                                                                       
  """                                                                           
  Rscript ${params.script_dir}/assign_observed_seqs_to_ref_seqs.R ${params.script_dir} $tidy_table $blast_output $metadata ${params.targets}
  """             
}


/* 
 This process aligns all of the unique sequences reported by dada2
 to the set of expected reference sequences using blastn.
*/
process blast_unassigned_sequences {
  publishDir "${params.outdir}", mode: 'link'

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

process assign_non_ref_seqs {
  publishDir "${params.outdir}", mode: 'link'

  input:
  path(blast_output) from post_blast_unassigned_ch
  path(unassigned_sequences) from post_blast_unassigned_seq_ch

  output:
  // path("unassigned_sequences.fasta") into post_assign_to_refseq_ch
  // path("identified_targets.xlsx")
  // path("identified_targets.tsv")

  script:                                                                       
  """                                                                           
  Rscript ${params.script_dir}/assign_non_ref_seqs.R ${params.script_dir} $unassigned_sequences $blast_output
  """             
}

