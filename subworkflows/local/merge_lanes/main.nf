/*
================================================================================
                                MERGING
================================================================================
*/

include { SAMBAMBA_MERGE }              from '../../../modules/nf-core/sambamba/merge/main'

workflow MERGE_LANES {
    take:
        take: bam_bwa

    main:
    ch_versions = Channel.empty()
    
        bam_bwa
            .map{meta, bam, bai ->
                meta.id = meta.patient + "_" + meta.sample
                [[meta.patient, meta.sample, meta.id, meta.gender, meta.status], bam, bai] }
            .groupTuple()
            .branch{
                single:   it[1].size() == 1
                multiple: it[1].size() > 1
            }.set{ bam_bwa_to_sort }
        
        bam_multiple = bam_bwa_to_sort.multiple
                                .map{
                                    info,bam,bai ->
                                    def meta = [:]
                                    meta.patient = info[0]
                                    meta.sample  = info[1]
                                    meta.id      = info[2]
                                    meta.gender  = info[3]
                                    meta.status  = info[4]
                                    [meta,bam,bai]
                                }
        
        bam_single = bam_bwa_to_sort.single
                                .map{
                                    info,bam,bai ->
                                    def meta = [:]
                                    meta.patient = info[0]
                                    meta.sample  = info[1]
                                    meta.id      = info[2]
                                    meta.gender  = info[3]
                                    meta.status  = info[4]
                                    [meta,bam,bai]
                                }
        // STEP 1.5: MERGING AND INDEXING BAM FROM MULTIPLE LANES
        //bam_multiple.view()
        SAMBAMBA_MERGE(bam_multiple)
        bam          = bam_single.mix(SAMBAMBA_MERGE.out.bam)
        ch_versions  = ch_versions.mix(SAMBAMBA_MERGE.out.versions.first())
    emit:
        bam
        ch_versions
}
