include { paramsSummaryMultiqc   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../../subworkflows/local/utils_nfcore_lp_wgs_pipeline'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { MULTIQC                } from '../../../modules/nf-core/multiqc/main'

workflow REPORTING_MULTIQC {
    take:
    versions
    reports
    outdir
    multiqc_config
    multiqc_logo
    multiqc_methods_description

    main:
    ch_multiqc_files = Channel.empty()

    ch_version_yaml = softwareVersionsToYAML(versions)
        .collectFile(storeDir: "${outdir}/pipeline_info", name: 'lp_wgs_software_mqc_versions.yml', sort: true, newLine: true)

    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = multiqc_config ? Channel.fromPath(multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo          = multiqc_logo ? Channel.fromPath(multiqc_logo, checkIfExists: true) : Channel.empty()
    multiqc_summary_params   = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(multiqc_summary_params))
    multiqc_methods_template = multiqc_methods_description ? file(multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description   = Channel.value(methodsDescriptionText(multiqc_methods_template))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'lp_wgs_workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_version_yaml)
    ch_multiqc_files = ch_multiqc_files.mix(reports)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'lp_wgs_methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    report = MULTIQC.out.report.toList()
}
