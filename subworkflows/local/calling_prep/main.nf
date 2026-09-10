include { HMMCOPY_GCCOUNTER } from '../../../modules/nf-core/hmmcopy/gccounter/main'
include { HMMCOPY_READCOUNTER } from '../../../modules/local/hmmcopy/readcounter/main'
include { SAMTOOLS_VIEW } from '../../../modules/local/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_NVIEW } from '../../../modules/local/samtools/view/main'
include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'

workflow CALLING_PREP {
    take:
    ch_input_sample
    ch_mapped_bam
    fasta
    fasta_fai
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
            ch_analysis_input = SAMTOOLS_NVIEW.out.bam
            versions = versions.mix(SAMTOOLS_NVIEW.out.versions.first())
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
    
    picard_input_bam = ch_analysis_input
                                .map{ meta, bam , bai -> tuple(meta, bam)}
    
    picard_fasta = fasta.join(fasta_fai).first()
    PICARD_MARKDUPLICATES ( picard_input_bam, picard_fasta)

    dups_marked = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai)
    HMMCOPY_READCOUNTER( dups_marked, fasta )

    emit:
    analysis_input = dups_marked
    gc_wig = ch_gc_wig
    readcounter_wig = HMMCOPY_READCOUNTER.out.wig
    versions
}
