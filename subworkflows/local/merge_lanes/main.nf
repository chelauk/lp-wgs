/*
================================================================================
                                MERGING
================================================================================
*/

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { SAMBAMBA_MERGE } from '../../../modules/local/sambamba/merge/main.nf'

workflow MERGE_LANES {
    take:
    ch_bam_bwa

    main:
    ch_bam_grouped = ch_bam_bwa
        .map { meta, bam ->
            def sampleMeta = meta.subMap(['patient', 'sample', 'status']) + [id: "${meta.patient}_${meta.sample}"]
            [sampleMeta, bam]
        }
        .groupTuple()
        .branch { _meta, bams ->
            single: bams.size() == 1
            multiple: true
        }

    ch_bam_single = ch_bam_grouped.single.map { meta, bams -> [meta, bams[0]] }
    // STEP 1.5: MERGING AND INDEXING BAM
    SAMBAMBA_MERGE(ch_bam_grouped.multiple) // already [meta, [bam1, bam2, ...]]
    ch_bam_unindexed = ch_bam_single.mix(SAMBAMBA_MERGE.out.bam)
    SAMTOOLS_INDEX(ch_bam_unindexed)
    
    emit:
    ch_bam_unindexed.join(SAMTOOLS_INDEX.out.index)
}
