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
 These fastq files represent the main input to this workflow
 
 Expecting files with _R1 or _R2 in their names corresponding to paired-end reads
*/


Channel
    .fromFilePairs("${params.fastq_dir}/${params.fastq_pattern}",
                   size: 2,
                   // checkIfExists: true,
                   maxDepth: 1)
    // this map step first strips off any _L### text from the end of the sample ID
    // and then any _S## text from the end
    // this is text added onto IDs from samplesheet that are added by Illumina bcl2fastq
    .map { untrimmed_sample_id, fastq  ->
           def sample_id = untrimmed_sample_id.replaceFirst( /_L\d{3}$/, "")
           sample_id = sample_id.replaceFirst( /_S\d+$/, "")
           [sample_id, fastq] }
    .set {reads_ch}


/*
  collect sample IDs in a file
*/
reads_ch.map{sample_id, fastq -> sample_id}
  .collectFile(name: 'sample_ids.txt', newLine: true)
  .set{ sample_ids_file_ch }


/*
   this combinatorially mixes the fastq file pairs with the primer sequences
   to be trimmed, creating a new channel with all possible combinations of
   input fastq and primer pair
*/                 
    
primers_and_samples_ch = primers_ch.combine(reads_ch)


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

                                                                                
//                                                                              
include { GENERATE_REFSEQ_FASTA       } from '../modules/local/generate_refseq_fasta/main'
include { SETUP_INDEXES               } from '../modules/local/setup_indexes/main'
include { SETUP_PYTHON_ENVIRONMENT    } from '../modules/local/setup_python_env/main'
include { SETUP_R_DEPENDENCIES        } from '../modules/local/setup_R_dependencies/main'
include { VALIDATE_METADATA           } from '../modules/local/validate_metadata/main'


include { FASTQC                      } from '../modules/nf-core/fastqc/main'
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
    SETUP_R_DEPENDENCIES(R_install_pkg_script_ch, R_tar_dir_ch, params.R_lib_dir)
    ch_versions = ch_versions.mix ( SETUP_R_DEPENDENCIES.out.versions )      

    // Check for existence of metadata file and create channel 
    validate_metadata_script_ch = Channel.fromPath(params.validate_metadata_script, checkIfExists: true)
    metadata_ch = Channel.fromPath("${params.metadata}", checkIfExists: true)
    VALIDATE_METADATA(validate_metadata_script_ch, metadata_ch, sample_ids_file_ch)
    ch_versions = ch_versions.mix ( VALIDATE_METADATA.out.versions )      

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
                                                                                
    MULTIQC (                                                                   
        ch_multiqc_files.collect(),                                             
        ch_multiqc_config.toList(),                                             
        ch_multiqc_custom_config.toList(),                                      
        ch_multiqc_logo.toList()                                                
    )                                                                           
    multiqc_report = MULTIQC.out.report.toList()                                

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