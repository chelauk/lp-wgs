include { BWA_MEM                              } from '../../../modules/local/bwa/mem/main'
include { QC_TRIM                              } from '../../../subworkflows/local/qc_trim/main'
include { MERGE_LANES                          } from '../../../subworkflows/local/merge_lanes/main'
include { MOSDEPTH                             } from '../../../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTINSERTSIZEMETRICS      } from '../../../modules/nf-core/picard/collectinsertsizemetrics/main'
include { SAMTOOLS_VIEW                        } from '../../../modules/local/samtools/view/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../../modules/local/picard/collectalignmentsummarymetrics/main'
include { PUBLISH_MAPPED_BAM                   } from '../../../modules/local/publish_mapped_bam/main'
include { PICARD_MARKDUPLICATES                } from '../../../modules/nf-core/picard/markduplicates/main'

workflow MAPPING_QC {
    take:
    ch_input_sample
    bwa
    fasta
    fasta_fai
    dict
    chr_bed
    sort
    fastp_adapter_fasta
    filter_bam
    filter_bam_min
    filter_bam_max
    filter_status

    main:
    reports  = channel.empty()
    versions = channel.empty()

    QC_TRIM(ch_input_sample, fastp_adapter_fasta)
    versions = versions.mix(QC_TRIM.out.versions)
    reports  = reports.mix(QC_TRIM.out.reports)

    BWA_MEM(QC_TRIM.out.reads, bwa, fasta, sort)
    versions = versions.mix(BWA_MEM.out.versions_bwa)
    versions = versions.mix(BWA_MEM.out.versions_samtools)

    MERGE_LANES(BWA_MEM.out.bam)
    versions = versions.mix(MERGE_LANES.out.versions.first())

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
        versions = versions.mix(SAMTOOLS_VIEW.out.versions.first())
    }
    

    picard_input_bam = ch_mapped_bam
                                .map{ meta, bam , _bai -> tuple(meta, bam)}
    
    picard_fasta = fasta.join(fasta_fai).first()
    PICARD_MARKDUPLICATES ( picard_input_bam, picard_fasta)

    dups_marked = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai)

    PUBLISH_MAPPED_BAM(dups_marked)
    ch_mapped_bam = PUBLISH_MAPPED_BAM.out.bam

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(
        ch_mapped_bam,
        fasta,
        dict,
        filter_status
    )
    versions = versions.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.versions.first())
    reports  = reports.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics.collect { _meta, report -> report })

    ch_insert_metrics_input = ch_mapped_bam.map { meta, bam, _bai -> [meta, bam] }
    PICARD_COLLECTINSERTSIZEMETRICS(ch_insert_metrics_input)
    versions = versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions_picard)
    reports  = reports.mix(
        PICARD_COLLECTINSERTSIZEMETRICS.out.metrics.collect { _meta, report -> report },
        PICARD_COLLECTINSERTSIZEMETRICS.out.histogram.collect { _meta, report -> report }
    )

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
    versions = versions.mix(MOSDEPTH.out.versions_mosdepth)
    versions = versions.mix(MOSDEPTH.out.versions_gzip)
    reports  = reports.mix(
        MOSDEPTH.out.global_txt.collect { _meta, report -> report },
        MOSDEPTH.out.summary_txt.collect { _meta, report -> report },
        MOSDEPTH.out.regions_txt.collect { _meta, report -> report }
    )

    emit:
    bam = ch_mapped_bam
    reports
    versions
}
