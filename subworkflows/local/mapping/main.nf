include { BWA_MEM             } from '../../../modules/local/bwa/mem/main'
include { BAM_QC              } from '../../../subworkflows/local/bam_qc/main'
include { QC_TRIM             } from '../../../subworkflows/local/qc_trim/main'
include { MERGE_LANES         } from '../../../subworkflows/local/merge_lanes/main'
include { SAMTOOLS_VIEW       } from '../../../modules/local/samtools/view/main'
include { PUBLISH_MAPPED_BAM  } from '../../../modules/local/publish_mapped_bam/main'

workflow MAPPING {
    take:
    ch_input_sample
    bwa
    fasta
    fasta_fai
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
    
    BAM_QC(
        ch_mapped_bam,
        fasta,
        fasta_fai,
        chr_bed,
        filter_status
        )

    PUBLISH_MAPPED_BAM(BAM_QC.out.dups_marked)
    
    reports = QC_TRIM.out.multiqc_files
    
    emit:
    bam = BAM_QC.out.dups_marked
    multiqc_files = reports
}
