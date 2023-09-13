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
    params.fasta_fai,
    params.centromere,
    params.map_wig,
    params.dict
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

include { BWA_MEM                     } from '../modules/nf-core/bwa/mem/main'
include { QC_TRIM                     } from '../subworkflows/local/qc_trim/main'
include { MERGE_LANES                 } from '../subworkflows/local/merge_lanes/main'
include { MOSDEPTH                    } from '../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTINSERTSIZEMETRICS } from '../modules/nf-core/picard/collectinsertsizemetrics/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { HMMCOPY_GCCOUNTER           } from '../modules/nf-core/hmmcopy/gccounter/main'
include { HMMCOPY_READCOUNTER         } from '../modules/nf-core/hmmcopy/readcounter/main'
include { ICHORCNA_RUN                } from '../modules/nf-core/ichorcna/run/main'
include { SAMTOOLS_VIEW               } from '../modules/nf-core/samtools/view/main'
include { ACE                         } from '../modules/local/ace'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../modules/local/picard/collectalignmentmummarymetrics/main'

//
// value channels
//

map_bin = params.map_bin

//
// gather prebuilt indices
//
    dict                   = params.dict               ? Channel.fromPath(params.dict).collect()      : Channel.empty()
    fasta                  = params.fasta              ? Channel.fromPath(params.fasta).collect()     : Channel.empty()
    fasta_fai              = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect() : Channel.empty()
    bwa                    = params.bwa                ? Channel.fromPath(params.bwa).collect()       : Channel.empty()
    chr_bed                = params.chr_bed            ? Channel.fromPath(params.chr_bed).collect()       : Channel.empty()
    centromere             = params.centromere         ? Channel.fromPath(params.centromere).collect(): Channel.empty()
    if ( params.map_bin == '10kb' ) {
        gc_wig                = params.map_wig            ? Channel.fromPath("${params.map_wig}/gc_hg38_10kb.wig").collect()   : Channel.empty()
        map_wig               = params.map_wig            ? Channel.fromPath("${params.map_wig}/map_hg38_10kb.wig").collect()   : Channel.empty()
        pon_rds               = params.normal             ? Channel.fromPath("${params.map_wig}/HD_ULP_PoN_hg38_10kb_median_normAutosome_median.rds").collect() : Channel.value([]) // optional
    } else if ( params.map_bin == '50kb' ) {
        gc_wig                = params.map_wig            ? Channel.fromPath("${params.map_wig}/gc_hg38_50kb.wig").collect()   : Channel.empty()
        map_wig               = params.map_wig            ? Channel.fromPath("${params.map_wig}/map_hg38_50kb.wig").collect()   : Channel.empty()
        pon_rds               = params.normal             ? Channel.fromPath("${params.map_wig}/HD_ULP_PoN_hg38_50kb_median_normAutosome_median.rds").collect() : Channel.value([]) // optional Channel.empty()
    } else if ( params.map_bin == '500kb') {
        gc_wig                = params.map_wig            ? Channel.fromPath("${params.map_wig}/gc_hg38_500kb.wig").collect()   : Channel.empty()
        map_wig               = params.map_wig            ? Channel.fromPath("${params.map_wig}/map_hg38_500kb.wig").collect()   : Channel.empty()
        pon_rds               = params.normal             ? Channel.fromPath("${params.map_wig}/HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds").collect() : Channel.value([]) // optional Channel.empty()
    } else  if ( params.map_bin == '1000kb' ){
        gc_wig                = params.map_wig            ? Channel.fromPath("${params.map_wig}/gc_hg38_1000kb.wig").collect()   : Channel.empty()
        map_wig               = params.map_wig            ? Channel.fromPath("${params.map_wig}/map_hg38_1000kb.wig").collect()   : Channel.empty()
        pon_rds               = params.normal             ? Channel.fromPath("${params.map_wig}/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds").collect() : Channel.value([]) // optional Channel.empty()
    }

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
    mosdepth_reports = Channel.empty()

    //
    // FASTQC, FASTP and BWA
    //

    sort_bam = true
    // Gather index for mapping given the chosen aligner
    ch_map_index = bwa
    if ( params.step == 'fastq' ) {
        // Create input channel
        fastq_input = ch_input_sample
        QC_TRIM ( fastq_input, ch_map_index, sort_bam)
        ch_versions = ch_versions.mix(QC_TRIM.out.ch_versions)
        ch_reports  = ch_reports.mix(QC_TRIM.out.ch_reports)
        BWA_MEM( QC_TRIM.out.reads,   ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort_bam) // If aligner is bwa-mem
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
        MERGE_LANES ( BWA_MEM.out.bam )
        ch_versions = ch_versions.mix(MERGE_LANES.out.ch_versions.first())
        ch_bam_input = MERGE_LANES.out.bam
    } else if ( params.step == 'bam'  &&  params.filter_bam == null ){
        ch_bam_input = ch_input_sample
    } else if ( params.step == 'bam'  &&  params.filter_bam != null ){
        ch_filter_input = ch_input_sample
        SAMTOOLS_VIEW ( ch_filter_input )
        ch_bam_input = SAMTOOLS_VIEW.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
    }

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS ( ch_bam_input , fasta, dict)
    ch_versions = ch_versions.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.versions.first())
    ch_reports  = ch_reports.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics.collect{meta, report -> report})
    PICARD_COLLECTINSERTSIZEMETRICS ( ch_bam_input )
    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())
    ch_reports  = ch_reports.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.size_metrics.collect{meta, report -> report})

    MOSDEPTH(
        ch_bam_input,
        chr_bed,
        fasta.map{ it -> [[id:it[0].baseName], it] }
        )

    mosdepth_reports = mosdepth_reports.mix(MOSDEPTH.out.global_txt,
                                            MOSDEPTH.out.regions_txt)
    ch_reports  = ch_reports.mix(mosdepth_reports.collect{meta, report -> report})
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    // run hmmcopygccounter
    if ( params.call_gc ) {
    HMMCOPY_GCCOUNTER(
        fasta.map{ it -> [[id:it[0].baseName], it] },
        map_bin
        )
    gc_wig = HMMCOPY_GCCOUNTER.out.wig
    ch_versions = ch_versions.mix(HMMCOPY_GCCOUNTER.out.versions)
    } else {
        gc_wig = gc_wig
    }

    // run hmmcopyreadcounter
    HMMCOPY_READCOUNTER(
        ch_bam_input
        )
    ch_versions = ch_versions.mix(HMMCOPY_READCOUNTER.out.versions)


    panel_of_normals = pon_rds
    normal_wig = []
    // run ichorcna
    ICHORCNA_RUN(
        HMMCOPY_READCOUNTER.out.wig,
        normal_wig,
        gc_wig,
        map_wig,
        panel_of_normals,
        centromere
        )
    ch_versions= ch_versions.mix(ICHORCNA_RUN.out.versions)

    // run ACE
    ACE(ch_bam_input)
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
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect())
    
    ch_multiqc_files = ch_multiqc_files.collect()
    
    MULTIQC (
        ch_multiqc_files
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

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
    // check file is not empty
    file(csv_file).withReader('UTF-8') { reader ->
        if (reader.readLine() == null) {
            log.error "CSV file is empty"
            return null
        }
    }
    // read csv file
    Channel.of(csv_file).splitCsv(header: true)
    .map { row ->
        if (!(row.patient && row.sample)) log.warn "Missing or unknown field in csv file header"
        [[row.patient.toString(), row.sample.toString()], row]
    }
    .groupTuple()
    .map { meta, rows ->
        size = rows.size()
        [rows, size]
        }.transpose()
          //A Transpose Function takes a collection of columns and returns a collection of rows.
          //The first row consists of the first element from each column. Successive rows are constructed similarly.
          //def result = [['a', 'b'], [1, 2], [3, 4]].transpose()
          //assert result == [['a', 1, 3], ['b', 2, 4]]
            .map{
            row, num_lanes ->
            def meta = [:]
            if (row.patient) meta.patient = row.patient.toString()
            if (row.sample)  meta.sample  = row.sample.toString()
            if (row.sex)  meta.sex  = row.sex.toString()
                else meta.sex = "NA"
            if (row.status)  meta.status  = row.status.toString()
                else meta.status = 0
            if (row.fastq_1) {
                meta.patient  = row.patient.toString()
                meta.sample   = row.sample.toString()
                if (row.lane) {
                    meta.lane = row.lane.toString()
                    meta.id   = "${row.patient}_${row.sample}_${row.lane}"
                } else {
                    meta.id = "${row.patient}_${row.sample}"
                    }
                def fastq_1 = file(row.fastq_1, checkIfExists: true)
                def fastq_2 = file(row.fastq_2, checkIfExists: true)
                def CN      = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
                def flowcell    = flowcellLaneFromFastq(fastq_1)
                //Don't use a random element for ID, it breaks resuming
                def read_group  = "\"@RG\\tID:${flowcell}.${row.sample}.${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
                meta.read_group  = read_group.toString()
                meta.data_type   = "fastq"
                if (params.step == 'fastq') return [meta, [fastq_1, fastq_2]]
                else {
                    log.error "Samplesheet contains fastq files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
                    System.exit(1)
                    }
                } else if (row.bam) {
                    meta.id = "${row.patient}_${row.sample}"
                    def bam = file(row.bam, checkIfExists: true)
                    def bai = file(row.bai, checkIfExists: true)
                    meta.data_type  = 'bam'
                    if (!(params.step == 'fastq')) return [meta, bam, bai]
                    else {
                        log.error "Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
                        System.exit(1)
                    }
                }
            }
    }

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
        // Check if running in stub mode
    if (workflow.stubRun) {
        // Return default value when in -stub mode
        return "DEFAULT_FLOWCELL_LANE"
    }
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
