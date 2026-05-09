import WorkflowWgs

include { paramsSummaryMultiqc   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
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
    multiqc_summary_params   = WorkflowWgs.paramsSummaryMap(workflow, params, "nextflow_schema.json")
    ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(multiqc_summary_params))
    multiqc_methods_template = multiqc_methods_description ? file(multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description   = Channel.value(WorkflowWgs.methodsDescriptionText(workflow, multiqc_methods_template))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'lp_wgs_workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_version_yaml)
    ch_multiqc_files = ch_multiqc_files.mix(reports)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'lp_wgs_methods_description_mqc.yaml', sort: true))

    ch_multiqc_input = ch_multiqc_files
        .collect()
        .map { files -> [files: files] }
        .combine(ch_multiqc_config.toList().map { config -> [base_config: config] })
        .combine(ch_multiqc_custom_config.toList().map { config -> [custom_config: config] })
        .combine(ch_multiqc_logo.toList().map { logo -> [logo: logo] })
        .map { multiqc_input ->
            def input = multiqc_input.inject([:]) { acc, value -> acc + value }
            def files = input.files
            def base_config = input.base_config ?: []
            def custom_config = input.custom_config ?: []
            def logo = input.logo ?: []
            def configs = base_config + custom_config
            def logo_file = logo ? logo[0] : []
            [[id: 'lp_wgs'], files, configs, logo_file, [], []]
        }

    MULTIQC(
        ch_multiqc_input
    )

    emit:
    report = MULTIQC.out.report.map { meta, report -> report }.toList()
}

//
// Convert legacy versions.yml files into a MultiQC-ready YAML string.
//
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [k.tokenize(':')[-1], v] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Convert nf-core tuple-style version entries into a MultiQC-ready YAML string.
//
def processVersionsFromTuple(version) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def process = version[0].toString().tokenize(':')[-1]
    def tool = version[1].toString()
    def tool_version = version[2].toString()
    return yaml.dumpAsMap([(process): [(tool): tool_version]]).trim()
}

//
// Convert version entries from either versions.yml files or nf-core tuple outputs.
//
def processVersionEntry(version) {
    if (version instanceof Collection && version.size() == 3) {
        return processVersionsFromTuple(version)
    }
    return processVersionsFromYAML(version)
}

//
// Get workflow version for MultiQC software versions.
//
def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${workflow.manifest.version}
        Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format.
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
                .unique()
                .map { version -> processVersionEntry(version) }
                .unique()
                .mix(Channel.of(workflowVersionToYAML()))
}
