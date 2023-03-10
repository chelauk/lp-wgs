//
// QC Trimming and Alignment
//

include { FASTQC } from '../../modules/nf-core/fastqc/main'

workflow QC_TRIM_ALIGN {
    take:
    reads               // channel: [ val(meta), [ reads ] ]

    main:
    FASTQC(reads)
    ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
    ch_versions = ch_versions.mix(FASTQ.out.versions.first())

    save_trimmed_fail = false
    save_merged = false
    FASTP(
        reads,
        [],  // default adapter sequences
        save_trimmed_fail,
        save_merged
        )
    ch_reports  = ch_reports.mix(
                            FASTP.out.json.collect{meta, json -> json},
                            FASTP.out.html.collect{meta, html -> html}
                            )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    BWA_MEM(FASTP.out.)

    emit:
    map
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
