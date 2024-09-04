/*
================================================================================
                                MERGING
================================================================================
*/

include { SAMTOOLS_INDEX              } from '../../../modules/nf-core/samtools/index/main'
include { SAMBAMBA_MERGE }              from '../../../modules/local/sambamba/merge/main.nf'

workflow MERGE_LANES {
    take:
        take: bam_bwa

    main:
    versions = Channel.empty()
    
        bam_bwa
            .map{meta, bam ->
                meta.id = meta.patient + "_" + meta.sample
                [[meta.patient, meta.sample, meta.id, meta.gender, meta.status], bam ] }
            .groupTuple()
            .branch{
                single:   it[1].size() == 1
                multiple: it[1].size() > 1
            }.set{ bam_bwa_to_sort }
        
        bam_multiple = bam_bwa_to_sort.multiple
                                .map{
                                    info,bam ->
                                    def meta = [:]
                                    meta.patient = info[0]
                                    meta.sample  = info[1]
                                    meta.id      = info[2]
                                    meta.gender  = info[3]
                                    meta.status  = info[4]
                                    [meta,bam]
                                }
        
        bam_single = bam_bwa_to_sort.single
                                .map{
                                    info,bam ->
                                    def meta = [:]
                                    meta.patient = info[0]
                                    meta.sample  = info[1]
                                    meta.id      = info[2]
                                    meta.gender  = info[3]
                                    meta.status  = info[4]
                                    [meta,bam]
                                }
        // STEP 1.5: MERGING AND INDEXING BAM FROM MULTIPLE LANES
        //bam_multiple.view()
        SAMTOOLS_INDEX(bam_single)
        bam_single = bam_single.join(SAMTOOLS_INDEX.out.bai)
        SAMBAMBA_MERGE(bam_multiple)
        bam          = bam_single.mix(SAMBAMBA_MERGE.out.bam)
        versions  = versions.mix(SAMBAMBA_MERGE.out.versions.first())
    emit:
        bam
        versions
}
