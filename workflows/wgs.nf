/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowWgs.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.bwa,
    params.fasta,
    params.fasta_fai
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
ch_input_sample = extract_csv(file(params.input, checkIfExists: true ))

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { QC_TRIM_ALIGN               } from '../subworkflows/local/qc_trim_align'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { ACE                         } from '../modules/local/ace'

//
// gather prebuilt indices
//
    fasta                  = params.fasta              ? Channel.fromPath(params.fasta).collect()     : Channel.empty()
    fasta_fai              = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect() : Channel.empty()
    bwa                    = params.bwa                ? Channel.fromPath(params.bwa).collect()       : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow WGS {
    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
    ch_versions = Channel.empty()
    ch_version_yaml = Channel.empty()

    //
    // FASTQC, FASTP and BWA
    //

    sort_bam = true
    // Gather index for mapping given the chosen aligner
    ch_map_index = bwa
    if ( params.step == 'mapping' ) {
        // Create input channel
        fastq_input = ch_input_sample
        QC_TRIM_ALIGN ( fastq_input, ch_map_index, sort_bam)
        ch_versions = ch_versions.mix(QC_TRIM_ALIGN.out.ch_versions)
        ch_reports  = ch_reports.mix(QC_TRIM_ALIGN.out.ch_reports)
        ch_ace_input = QC_TRIM_ALIGN.out.bam
    } else if ( params.step == 'ace' ) {
        ch_ace_input = ch_input_sample
    }

    // run ACE
    ACE(ch_ace_input)
    ch_versions = ch_versions.mix(ACE.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowWgs.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowWgs.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.patient && row.sample)) log.warn "Missing or unknown field in csv file header"
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing
        def meta = [:]

        //TODO since it is mandatory: error/warning if not present?
        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no gender specified, gender is not considered
        // gender is only mandatory for somatic CNV
        if (row.gender) meta.gender = row.gender.toString()
        else meta.gender = "NA"

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        // mapping with fastq
        if (row.fastq_1) {
            meta.patient    = row.patient.toString()
            meta.sample     = row.sample.toString()
            if (row.lane){
                meta.lane       = row.lane.toString()
                meta.id         = "${row.patient}_${row.sample}_${row.lane}"}
                else {
                    meta.id     = "${row.patient}_${row.sample}"
                }
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.patient}_${row.sample}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.patient}_${row.sample}\\tPL:ILLUMINA\""
            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'fastq'

            fqs = [] // list to fill with all arguments matching regex below
            map_fqs = row.findAll{k,v -> k.matches(~/^fastq(_[0-9]+)?/)} // find fqs in row
            map_fqs = map_fqs.sort{it -> (it.key.replaceAll(/fastq(_)?/, "") ?: 0).toInteger() } // sort matches
            for (fq in map_fqs) {
            fqs.add(file(fq.value, checkIfExists: true))
            }

            return [meta, fqs]

        // start from BAM
        } else if (row.patient && row.bam) {
            meta.patient    = row.patient.toString()
            meta.sample     = row.sample.toString()
            meta.id         = "${row.patient}_${row.sample}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def read_group  = "\"@RG\\tID:${row.patient}_${row.sample}\\tLB:${row.patient}_${row.sample}\\tSM:${row.patient}_${row.sample}\\tPL:ILLUMINA\""
            meta.read_group = read_group.toString()
            meta.data_type  = "bam"
            return [meta, bam]
        } else {
            log.warn "Missing2 or unknown field in csv file header"
        }
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
