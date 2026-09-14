include { BWA_MEM                               } from '../../../modules/local/bwa/mem/main'
include { QC_TRIM                               } from '../../../subworkflows/local/qc_trim/main'
include { MERGE_LANES                           } from '../../../subworkflows/local/merge_lanes/main'
include { MOSDEPTH                              } from '../../../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../../modules/nf-core/picard/collectalignmentsummarymetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS       } from '../../../modules/nf-core/picard/collectinsertsizemetrics/main'
include { PICARD_MARKDUPLICATES                 } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_VIEW                         } from '../../../modules/local/samtools/view/main'
//include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../../modules/local/picard/collectalignmentsummarymetrics/main'
include { PUBLISH_MAPPED_BAM                    } from '../../../modules/local/publish_mapped_bam/main'

workflow MAPPING_QC {
    take:
    ch_input_sample
    bwa
    fasta
    fasta_fai
//    dict
    chr_bed
    sort
    fastp_adapter_fasta
    filter_bam
    filter_bam_min
    filter_bam_max
    filter_status

    main:
    QC_TRIM(ch_input_sample, fastp_adapter_fasta)

    BWA_MEM(QC_TRIM.out.reads, bwa, fasta, sort)

    MERGE_LANES(BWA_MEM.out.bam)

    if (!filter_bam) {
        ch_mapped_bam = MERGE_LANES.out.bam.map { meta, bam, bai ->
            [meta + [filter_status: filter_status], bam, bai]
        }
    } else {
        ch_filter_input = MERGE_LANES.out.bam.map { meta, bam, bai ->
            [meta + [filter_status: filter_status], bam, bai]
        }
        SAMTOOLS_VIEW(ch_filter_input, filter_bam_min, filter_bam_max)
        ch_mapped_bam = SAMTOOLS_VIEW.out.bam
    }
    

    picard_input_bam = ch_mapped_bam
                                .map{ meta, bam , _bai -> tuple(meta, bam)}
    // We use .first() to convert the single FASTA/FAI result from a queue channel into a reusable value channel 
    picard_fasta = fasta.join(fasta_fai).first()
    PICARD_MARKDUPLICATES ( picard_input_bam, picard_fasta)

    dups_marked = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai)

    PUBLISH_MAPPED_BAM(dups_marked)
    ch_mapped_bam = PUBLISH_MAPPED_BAM.out.bam

    // COLLECTALIGNMENTSUMMARYMETRICS requires [meta, fasta], without the FAI.
    alignment_summary_fasta = picard_fasta.map {
                                meta, fasta_file, _fai -> tuple(meta, fasta_file)
                               }

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(
        PICARD_MARKDUPLICATES.out.bam,
        alignment_summary_fasta
    )

    PICARD_COLLECTINSERTSIZEMETRICS(PICARD_MARKDUPLICATES.out.bam)

    ch_mosdepth_input = ch_mapped_bam
        .combine(chr_bed)
        .map { meta, bam, bai, bed ->
            def interval_bed = bed instanceof List ? bed[0] : bed
            [meta + [filter_status: filter_status], bam, bai, interval_bed]
        }

    MOSDEPTH(
        ch_mosdepth_input,
        fasta,
        []
    )
    reports = QC_TRIM.out.multiqc_files
    reports  = reports
                  .mix(PICARD_MARKDUPLICATES.out.metrics.collect { _meta, metrics -> metrics })
                  .mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics.collect { _meta, report -> report })
                  .mix(PICARD_COLLECTINSERTSIZEMETRICS.out.metrics.collect { _meta, report -> report })
                  .mix(
                       MOSDEPTH.out.global_txt.collect { _meta, report -> report },
                       MOSDEPTH.out.summary_txt.collect { _meta, report -> report },
                       MOSDEPTH.out.regions_txt.collect { _meta, report -> report })

    emit:
    bam = ch_mapped_bam
    multiqc_files = reports
}
