//
// Subworkflow that uses the nf-schema plugin to render help text and parameter summary
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-VALIDATION PLUGIN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsHelp         } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTILS_NFSCHEMA_PLUGIN {
    take:
    input_workflow
    validate_params
    schema_filename
    help
    help_full
    show_hidden
    pre_help_text
    post_help_text

    main:
    
    // Default values for optional strings
    pre_help_text  = pre_help_text  ?: ''
    post_help_text = post_help_text ?: ''
    
    log.debug("Using schema file: ${schema_filename}")
    
    //
    // Print help and exit before validating required parameters
    //
    if (help || help_full) {
        help_options = [
            beforeText : pre_help_text,
            afterText  : post_help_text,
            command    : "nextflow run ${workflow.manifest.name} --input samplesheet.csv --outdir results -profile <docker/singularity>",
            showHidden : show_hidden,
            fullHelp   : help_full
        ]
    
        if (schema_filename) {
            help_options.parameters_schema = schema_filename
        }
    
        log.info paramsHelp(
            help_options,
            (help instanceof String && help != 'true') ? help : ''
        )
    
        exit 0
    }
    
    //
    // Print parameters that differ from their schema defaults
    //
    summary_options = [:]
    
    if (schema_filename) {
        summary_options.parameters_schema = schema_filename
    }
    
    log.info pre_help_text
    log.info paramsSummaryLog(summary_options, input_workflow)
    log.info post_help_text
    
    //
    // Validate parameters against the JSON schema
    //
    if (validate_params) {
        validate_options = [:]
    
        if (schema_filename) {
            validate_options.parameters_schema = schema_filename
        }
    
        validateParameters(validate_options)
    }   
 
    emit:
    dummy_emit = true
}
