process PUBLISH_MAPPED_BAM {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path("*.bam"), path("*.{bai,csi,crai}"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    true
    """
}
