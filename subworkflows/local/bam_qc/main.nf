include { MOSDEPTH                              } from '../../../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../../modules/nf-core/picard/collectalignmentsummarymetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS       } from '../../../modules/nf-core/picard/collectinsertsizemetrics/main'
include { PICARD_MARKDUPLICATES                 } from '../../../modules/nf-core/picard/markduplicates/main'


workflow BAM_QC {

    take:
    ch_mapped_bam
    fasta
    fasta_fai
    chr_bed
    filter_status
    
    main:
    reports = channel.empty()
    picard_input_bam = ch_mapped_bam
                                .map{ meta, bam , _bai -> tuple(meta, bam)}
    // We use .first() to convert the single FASTA/FAI result from a queue channel into a reusable value channel 
    picard_fasta = fasta.join(fasta_fai).first()
    PICARD_MARKDUPLICATES ( picard_input_bam, picard_fasta)
    dups_marked = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai)
    
    // COLLECTALIGNMENTSUMMARYMETRICS requires [meta, fasta], without the FAI.
    alignment_summary_fasta = picard_fasta.map {
                                meta, fasta_file, _fai -> tuple(meta, fasta_file)
                               }

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(
        PICARD_MARKDUPLICATES.out.bam,
        alignment_summary_fasta
    )

    PICARD_COLLECTINSERTSIZEMETRICS(PICARD_MARKDUPLICATES.out.bam)

    ch_mosdepth_input = dups_marked
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
    reports  = reports
                  .mix(PICARD_MARKDUPLICATES.out.metrics.collect { _meta, metrics -> metrics })
                  .mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics.collect { _meta, report -> report })
                  .mix(PICARD_COLLECTINSERTSIZEMETRICS.out.metrics.collect { _meta, report -> report })
                  .mix(
                       MOSDEPTH.out.global_txt.collect { _meta, report -> report },
                       MOSDEPTH.out.summary_txt.collect { _meta, report -> report },
                       MOSDEPTH.out.regions_txt.collect { _meta, report -> report })

    emit:
    dups_marked
    multiqc_files = reports

}
