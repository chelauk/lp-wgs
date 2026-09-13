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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTILS_NFSCHEMA_PLUGIN {
    take:
    pre_help_text    //  string: string to be printed before help text and summary log
    post_help_text   //  string: string to be printed after help text and summary log
    validate_params  // boolean: validate parameters
    schema_filename  //    path: JSON schema file, null to use default value

    main:

    log.debug("Using schema file: ${schema_filename}")

    // Default values for strings
    pre_help_text    = pre_help_text    ?: ''
    post_help_text   = post_help_text   ?: ''

    //
    // Print parameter summary to stdout
    //
    log.info(pre_help_text + paramsSummaryLog(workflow, parameters_schema: schema_filename) + post_help_text)

    //
    // Validate parameters relative to the parameter JSON schema
    //
    if (validate_params) {
        validateParameters(parameters_schema: schema_filename)
    }

    emit:
    true
}
