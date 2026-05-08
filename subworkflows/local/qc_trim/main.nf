//
// QC Trimming and Alignment
//

include { FASTQC  } from '../../../modules/nf-core/fastqc/main'
include { FASTP   } from '../../../modules/nf-core/fastp/main'

workflow QC_TRIM {
    take:
    ch_reads // channel: [ val(meta), fastq_1, fastq_2 ]
    fastp_adapter_fasta

    main:
    versions = Channel.empty()
    reports  = Channel.empty()

    ch_reads
        .map { meta, fastq_1, fastq_2 -> [meta, [fastq_1, fastq_2]] }
        .set { ch_reads_for_qc }

    FASTQC(ch_reads_for_qc)
    reports  = reports.mix(FASTQC.out.zip.collect { meta, logs -> logs })
    versions = versions.mix(FASTQC.out.versions.first())

    adapter_fasta = fastp_adapter_fasta ? file(fastp_adapter_fasta, checkIfExists: true) : []

    ch_fastp_input = ch_reads_for_qc.map { meta, reads ->
        tuple(meta, reads, adapter_fasta)
    }

    discard_trimmed_pass = false
    save_trimmed_fail = false
    save_merged = false

    FASTP(
        ch_fastp_input,
        discard_trimmed_pass,
        save_trimmed_fail,
        save_merged
    )
    reports  = reports.mix(
        FASTP.out.json.collect { meta, json -> json },
        FASTP.out.html.collect { meta, html -> html }
    )

    emit:
    reads = FASTP.out.reads // channel: [ val(meta), [ reads ] ]
    versions    // channel: fastqc, fastp versions
    reports    // channel: fastqc, fastp reports
}
