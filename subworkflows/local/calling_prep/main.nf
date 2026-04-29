include { HMMCOPY_GCCOUNTER } from '../../../modules/nf-core/hmmcopy/gccounter/main'
include { HMMCOPY_READCOUNTER } from '../../../modules/nf-core/hmmcopy/readcounter/main'
include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/'
include { SAMTOOLS_NVIEW } from '../../../modules/nf-core/samtools/view_nanopore/main'

workflow CALLING_PREP {
    take:
    ch_input_sample
    ch_mapped_bam
    fasta
    gc_wig
    step
    tech
    filter_bam
    filter_bam_min
    filter_bam_max
    call_gc
    bin_size

    main:
    versions = Channel.empty()

    ch_analysis_input = step == 'calling' ? ch_input_sample : ch_mapped_bam

    if (tech == 'illumina') {
        if (step == 'calling' && filter_bam) {
            ch_filter_input = ch_analysis_input
            SAMTOOLS_VIEW(ch_filter_input, filter_bam_min, filter_bam_max)
            ch_analysis_input = SAMTOOLS_VIEW.out.bam
            versions = versions.mix(SAMTOOLS_VIEW.out.versions.first())
        }
    } else if (tech == 'nanopore') {
        if (step == 'calling' && filter_bam) {
            ch_filter_input = ch_analysis_input
            SAMTOOLS_NVIEW(ch_filter_input, filter_bam_min, filter_bam_max)
            ch_analysis_input = SAMTOOLS_NVIEW.out.bam.map { meta, bam, bai -> [meta, [bam, bai]] }
            versions = versions.mix(SAMTOOLS_NVIEW.out.versions.first())
        }
    } else {
        exit 1, "Unsupported sequencing technology '${tech}'. Expected one of: illumina, nanopore."
    }

    if (call_gc) {
        HMMCOPY_GCCOUNTER(fasta, bin_size)
        ch_gc_wig = HMMCOPY_GCCOUNTER.out.wig
        versions = versions.mix(HMMCOPY_GCCOUNTER.out.versions)
    } else {
        ch_gc_wig = gc_wig
    }

    HMMCOPY_READCOUNTER(ch_analysis_input)
    versions = versions.mix(HMMCOPY_READCOUNTER.out.versions)

    emit:
    analysis_input = ch_analysis_input
    gc_wig = ch_gc_wig
    readcounter_wig = HMMCOPY_READCOUNTER.out.wig
    versions
}
