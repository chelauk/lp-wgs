/*
================================================================================
                                MERGING
================================================================================
*/

include { SAMTOOLS_INDEX              } from '../../../modules/nf-core/samtools/index/main'
include { SAMBAMBA_MERGE }              from '../../../modules/local/sambamba/merge/main.nf'

workflow MERGE_LANES {
    take:
    ch_bam_bwa

    main:
    versions = Channel.empty()

    ch_bam_bwa
        .map { meta, bam ->
            def mergedMeta = meta + [id: "${meta.patient}_${meta.sample}"]
            [mergedMeta, bam]
        }
        .map { meta, bam -> [[meta.patient, meta.sample], [meta, bam]] }
        .groupTuple()
        .branch { sample_key, grouped_records ->
            single: grouped_records.size() == 1
            multiple: grouped_records.size() > 1
        }
        .set { ch_bam_grouped }

    ch_bam_single = ch_bam_grouped.single
        .map { sample_key, grouped_records -> grouped_records[0] }

    ch_bam_multiple = ch_bam_grouped.multiple
        .map { sample_key, grouped_records ->
            def meta = grouped_records[0][0]
            def bams = grouped_records.collect { record -> record[1] }
            [meta, bams]
        }

    // STEP 1.5: MERGING AND INDEXING BAM FROM MULTIPLE LANES
    SAMTOOLS_INDEX(ch_bam_single)
    ch_bam_single = ch_bam_single.join(SAMTOOLS_INDEX.out.bai)
    SAMBAMBA_MERGE(ch_bam_multiple)
    bam = ch_bam_single.mix(SAMBAMBA_MERGE.out.bam)
    versions = versions.mix(SAMBAMBA_MERGE.out.versions.first())

    emit:
    bam
    versions
}
