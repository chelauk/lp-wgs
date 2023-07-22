//
// QC Trimming and Alignment
//

include { FASTQC  } from '../../../modules/nf-core/fastqc/main'
include { FASTP   } from '../../../modules/nf-core/fastp/main'

workflow QC_TRIM {
    take:
    ch_reads               // channel: [ val(meta), [ reads ] ]
    ch_map_index           // channel: [mandatory] mapping index
    sort                   // boolean: [mandatory] true -> sort, false -> don't sort

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    FASTQC(ch_reads)
    ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    save_trimmed_fail = false
    save_merged = false
    FASTP(
        ch_reads,
        [],  // default adapter sequences
        save_trimmed_fail,
        save_merged
        )
    ch_reports  = ch_reports.mix(
                            FASTP.out.json.collect{meta, json -> json},
                            FASTP.out.html.collect{meta, html -> html}
                            )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    // BWA_MEM(FASTP.out.reads,   ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort) // If aligner is bwa-mem
    // Get the bam files from the aligner

    emit:
    reads          = FASTP.out.reads // channel: [ [meta], reads ]
    ch_versions    // channel: fastqc, fastp versions
    ch_reports    // channel: fastqc, fastp reports
}
