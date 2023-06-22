/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS                                                             
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)            

// This list includes a list of files or paths that are required
// to exist.  Check that they exist and fail if not.
checkPathParamList = [
  params.fastq_dir,
  params.script_dir,
  params.targets,
  params.metadata,
  params.primers
]
// log.info("Checking for required input paths and files...")
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
                                                                                
// Validate input parameters                                                    
WorkflowMPAS.initialise(params, log)                                           



/*
  Read in the targets tsv-format file that describes the expected target sequences.
  and create channel
*/
Channel
    .fromPath(params.targets, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { targets_ch }

/*
  A channel containing the path of the tsv file specifying surveillance targets
  
  This channel differs from targets_ch in that targets_ch contains a groovy
  data structure from the parsed file.  This ch contains the file's path.
*/
Channel
     .fromPath(params.targets, type: 'file', checkIfExists: true)
     .set{ targets_file_ch }

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
  A channel containing the tsv file specifying surveillance report columns
*/
Channel
     .fromPath(params.surveillance_columns, type: 'file', checkIfExists: true)
     .set{surveillance_columns_ch}


                                                                                

/*
 These fastq files represent the main input to this workflow
 
 Expecting files with _R1 or _R2 in their names corresponding to paired-end reads
*/


Channel
    .fromFilePairs("${params.fastq_dir}/${params.fastq_pattern}",
                   size: 2,
                   maxDepth: 1)
    .map { untrimmed_sample_id, fastq  ->
           // strip off any _L### text from the end of the sample ID
           // and then any _S## text from the end
           // this is text added onto IDs from samplesheet that are added by Illumina bcl2fastq
           def sample_id = untrimmed_sample_id.replaceFirst( /_L\d{3}$/, "")
           sample_id = sample_id.replaceFirst( /_S\d+$/, "")

           // define a new empty map named meta for each sample                   
           // and populate it with id and single_end values                       
           // for compatibility with nf-core module expected parameters           
           // reads are just the list of fastq                                    
           def meta        = [:]      

           // this map gets rid of any _S\d+ at the end of sample IDs but leaves fastq
           // names alone.  E.g. strip _S1 from the end of a sample ID..          
           // This is typically sample #s from Illumina basecalling.              
           // could cause an issue if sequenced the same sample with              
           // multiple barcodes so was repeated on a sample sheet.                
           meta.id         = sample_id
                                                                                  
           // this pipeline only designed to work with paired-end data
           meta.single_end =  false

           // this last statement in the map closure is the return value          
           [ meta, fastq ] }     
      
    .set {reads_ch}

/*
  collect sample IDs in a file
*/

reads_ch.map{meta, reads -> meta.id}
  .collectFile(name: 'sample_ids.txt', newLine: true)
  .set{ sample_ids_file_ch }

/*
   this combinatorially mixes the fastq file pairs with the primer sequences
   to be trimmed, creating a new channel with all possible combinations of
   input fastq and primer pair
*/                 
    
primers_and_reads_ch = primers_ch.combine(reads_ch)


/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES                                                                
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
                                                                                
/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS                                           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              

                                                                                
// local                                                                            
include { GENERATE_REFSEQ_FASTA       } from '../subworkflows/generate_refseq_fasta'
include { SETUP_INDEXES               } from '../modules/local/setup_indexes/main'
include { SETUP_PYTHON_ENVIRONMENT    } from '../modules/local/setup_python_env/main'
include { SETUP_R_DEPENDENCIES        } from '../modules/local/setup_R_dependencies/main'
include { VALIDATE_METADATA           } from '../modules/local/validate_metadata/main'
include { TRIM_PRIMER_SEQS            } from '../modules/local/trim_primer_seqs/main'
include { COLLECT_CUTADAPT_OUTPUT     } from '../modules/local/collect_cutadapt_output/main'
include { DADA2                       } from '../modules/local/dada2/main'
include { TIDY_DADA_OUTPUT            } from '../modules/local/dada2/main'
include { COMPARE_OBSERVED_SEQS       } from '../modules/local/assign_sequences/main'
include { ASSIGN_OBSERVED_SEQS        } from '../modules/local/assign_sequences/main'
include { GENERATE_TREES              } from '../subworkflows/generate_trees'
include { SETUP_BLAST_DB_AND_TAX      } from '../subworkflows/setup_blast_db_and_tax'
include { BLAST_UNASSIGNED_SEQUENCES  } from '../subworkflows/blast_unassigned_sequences'
include { PREPEND_OUTPUT_FILENAMES    } from '../modules/local/prepend_filenames/main'

// modules from NF-CORE 
include { FASTQC as FASTQC_PRE        } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST       } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW                                                           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
// Info required for completion email and summary                               
def multiqc_report = []                                                         
                                                                                
workflow MPAS_PIPELINE {                                                                
                                                                                
    ch_versions = Channel.empty()                                               
                                                                                
    CUSTOM_DUMPSOFTWAREVERSIONS (                                               
        ch_versions.unique().collectFile(name: 'collated_versions.yml')         
    )                     

    GENERATE_REFSEQ_FASTA(targets_ch)

    SETUP_INDEXES (GENERATE_REFSEQ_FASTA.out.fasta)
    ch_versions = ch_versions.mix ( SETUP_INDEXES.out.versions )      

    python_req_ch = Channel.fromPath(params.python_requirements, checkIfExists: true)
    SETUP_PYTHON_ENVIRONMENT(params.python_venv_path, python_req_ch)
    ch_versions = ch_versions.mix ( SETUP_PYTHON_ENVIRONMENT.out.versions )      

    R_install_pkg_script_ch = Channel.fromPath(params.R_install_pkg_script, checkIfExists: true)
    R_tar_dir_ch = Channel.fromPath(params.R_tar_dir, checkIfExists: true)
    SETUP_R_DEPENDENCIES(R_install_pkg_script_ch, R_tar_dir_ch)
    ch_versions = ch_versions.mix ( SETUP_R_DEPENDENCIES.out.versions )      

    // Check for existence of metadata file and create channel 
    validate_metadata_script_ch = Channel.fromPath(params.validate_metadata_script, checkIfExists: true)
    metadata_ch = Channel.fromPath("${params.metadata}", checkIfExists: true)
    VALIDATE_METADATA(validate_metadata_script_ch, metadata_ch, sample_ids_file_ch)
    ch_versions = ch_versions.mix ( VALIDATE_METADATA.out.versions )      

    // run fastqc on input reads
    FASTQC_PRE ( reads_ch )
    ch_versions = ch_versions.mix ( FASTQC_PRE.out.versions )      

    // identify and trim expected primer sequences
    // this happens once per sample x primer combination
    TRIM_PRIMER_SEQS ( primers_and_reads_ch )
    ch_versions = ch_versions.mix ( TRIM_PRIMER_SEQS.out.versions )      
 
    // concatenate individually-trimmed outputs and concatenate per sample
    COLLECT_CUTADAPT_OUTPUT ( TRIM_PRIMER_SEQS.out.trimmed_reads.groupTuple())
    ch_versions = ch_versions.mix ( COLLECT_CUTADAPT_OUTPUT.out.versions )      

    // run fastqc on trimmed reads
    FASTQC_POST ( COLLECT_CUTADAPT_OUTPUT.out.trimmed_reads )
    ch_versions = ch_versions.mix ( FASTQC_POST.out.versions )      

    // run dada2 on trimmed reads
    // DADA2 ( COLLECT_CUTADAPT_OUTPUT.out.trimmed_reads.collect() )
    ch_all_trimmed_reads = COLLECT_CUTADAPT_OUTPUT.out.trimmed_reads
                            .map { meta, reads -> [ reads ] }                                   
                            .flatten().collect()
    DADA2(ch_all_trimmed_reads)
    ch_versions = ch_versions.mix ( DADA2.out.versions )      

    TIDY_DADA_OUTPUT(DADA2.out.seqtab)
    ch_versions = ch_versions.mix ( TIDY_DADA_OUTPUT.out.versions )      

    COMPARE_OBSERVED_SEQS(TIDY_DADA_OUTPUT.out.observed_sequences,
                          TIDY_DADA_OUTPUT.out.sequence_abundance_table,
                          GENERATE_REFSEQ_FASTA.out.fasta,
                          SETUP_INDEXES.out.db)
    ch_versions = ch_versions.mix ( COMPARE_OBSERVED_SEQS.out.versions )      

    ASSIGN_OBSERVED_SEQS (VALIDATE_METADATA.out.validated_metadata,
                          COMPARE_OBSERVED_SEQS.out.sequence_abundance_table,
                          COMPARE_OBSERVED_SEQS.out.blast_output,
                          SETUP_R_DEPENDENCIES.out.R_lib_dir,
                          surveillance_columns_ch,
                          targets_file_ch)
    ch_versions = ch_versions.mix ( ASSIGN_OBSERVED_SEQS.out.versions )      

    GENERATE_TREES (ASSIGN_OBSERVED_SEQS.out.surveillance_report,
                    SETUP_PYTHON_ENVIRONMENT.out.venv_path,
                    targets_file_ch)      
    ch_versions = ch_versions.mix ( GENERATE_TREES.out.versions )      

    SETUP_BLAST_DB_AND_TAX()
    ch_versions = ch_versions.mix ( SETUP_BLAST_DB_AND_TAX.out.versions )      

    BLAST_UNASSIGNED_SEQUENCES(ASSIGN_OBSERVED_SEQS.out.unassigned_sequences,
                               SETUP_BLAST_DB_AND_TAX.out.blast_db_dir,
                               SETUP_BLAST_DB_AND_TAX.out.blast_tax_dir,
                               SETUP_R_DEPENDENCIES.out.R_lib_dir)

    ch_versions = ch_versions.mix ( BLAST_UNASSIGNED_SEQUENCES.out.versions )      



    //                                                                          
    // MODULE: MultiQC                                                          
    //                                                                          
    workflow_summary    = WorkflowMPAS.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)                       
                                                                                
    methods_description    = WorkflowMPAS.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)                 
                                                                                
    ch_multiqc_files = Channel.empty()                                          
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect{it[1]}.ifEmpty([]))
                                                                                
    MULTIQC (                                                                   
        ch_multiqc_files.collect(),                                             
        ch_multiqc_config.toList(),                                             
        ch_multiqc_custom_config.toList(),                                      
        ch_multiqc_logo.toList()                                                
    )                                                                           
    multiqc_report = MULTIQC.out.report.toList()                                

    // prepend main output file names
    ch_main_output_files = Channel.empty()

    ch_main_output_files = ch_main_output_files
                             .mix(multiqc_report, 
                             ASSIGN_OBSERVED_SEQS.out.surveillance_report,
                             ASSIGN_OBSERVED_SEQS.out.all_data_csv,
                             ASSIGN_OBSERVED_SEQS.out.txt,
                             ASSIGN_OBSERVED_SEQS.out.pdf,
                             BLAST_UNASSIGNED_SEQUENCES.out.report)
                             .flatten()
    PREPEND_OUTPUT_FILENAMES(ch_main_output_files)

}                                                                               
                                                                                
/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY                                                
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
workflow.onComplete {                                                           
/*
    TODO: fix this emailing functionality
    if (params.email || params.email_on_fail) {                                 
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }                                                                           
*/
    NfcoreTemplate.summary(workflow, params, log)                               
    if (params.hook_url) {                                                      
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }                                                                           
}                                                                               
                                                                                
/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END                                                                     
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/  
