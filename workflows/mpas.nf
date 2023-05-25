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
  Check for existence of metadata file 
  and create channel 
*/
Channel
    .fromPath("${params.metadata}",
                   checkIfExists: true)
    .set {post_metadata_check_ch}


/*
  Read in the targets tsv-format file that describes the expected target sequences.
  and create channel
*/
Channel
    .fromPath(params.targets, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { targets_ch }
                                                                                

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

    //                                                                          
    // MODULE: MultiQC                                                          
    //                                                                          
    workflow_summary    = WorkflowDbtlc.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)                       
                                                                                
    methods_description    = WorkflowDbtlc.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
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
