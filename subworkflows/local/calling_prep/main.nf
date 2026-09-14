include { HMMCOPY_GCCOUNTER } from '../../../modules/nf-core/hmmcopy/gccounter/main'
include { HMMCOPY_READCOUNTER } from '../../../modules/local/hmmcopy/readcounter/main'
include { SAMTOOLS_VIEW } from '../../../modules/local/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_NVIEW } from '../../../modules/local/samtools/view/main'
//include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'

workflow CALLING_PREP {
    take:
    ch_input_sample
    ch_mapped_bam
    fasta
//    fasta_fai
    gc_wig
    step
    tech
    filter_bam
    filter_bam_min
    filter_bam_max
    call_gc

    main:

    analysis_input = step == 'calling' ? ch_input_sample : ch_mapped_bam

    if (tech == 'illumina') {
        if (step == 'calling' && filter_bam) {
            ch_filter_input = analysis_input
            SAMTOOLS_VIEW(ch_filter_input, filter_bam_min, filter_bam_max)
            analysis_input = SAMTOOLS_VIEW.out.bam
        }
    } else if (tech == 'nanopore') {
        if (step == 'calling' && filter_bam) {
            ch_filter_input = analysis_input
            SAMTOOLS_NVIEW(ch_filter_input, filter_bam_min, filter_bam_max)
            analysis_input = SAMTOOLS_NVIEW.out.bam
        }
    } else {
        exit 1, "Unsupported sequencing technology '${tech}'. Expected one of: illumina, nanopore."
    }

    if (call_gc) {
        HMMCOPY_GCCOUNTER(fasta)
        ch_gc_wig = HMMCOPY_GCCOUNTER.out.wig
    } else {
        ch_gc_wig = gc_wig
    }
    
    HMMCOPY_READCOUNTER( analysis_input, fasta )

    emit:
    analysis_input
    gc_wig = ch_gc_wig
    readcounter_wig = HMMCOPY_READCOUNTER.out.wig
}
