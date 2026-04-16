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
    versions = Channel.empty()
    reports  = Channel.empty()
    
    ch_reads
        .map{ meta, fastq_1, fastq_2 -> [meta, [fastq_1,fastq_2]] }
        .view{"ch_reads before it goes into FASTP: $it"}
        .set{ ch_reads_for_qc }
    
    FASTQC(ch_reads_for_qc)
    reports  = reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
    versions = versions.mix(FASTQC.out.versions.first())

    save_trimmed_fail = false
    save_merged = false
	FASTP(
        ch_reads_for_qc,
        [],  // default adapter sequences
        save_trimmed_fail,
        save_merged
        )
    reports  = reports.mix(
                            FASTP.out.json.collect{meta, json -> json},
                            FASTP.out.html.collect{meta, html -> html}
                            )
    versions = versions.mix(FASTP.out.versions.first())
    // BWA_MEM(FASTP.out.reads,   ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort) // If aligner is bwa-mem
    // Get the bam files from the aligner

    //FASTP.out.reads.view{"ch_reads after it comes out of FASTP: $it"}
    emit:
    reads          = FASTP.out.reads // channel: [ [meta], reads ]
    versions    // channel: fastqc, fastp versions
    reports    // channel: fastqc, fastp reports
}
