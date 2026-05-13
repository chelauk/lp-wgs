process PUBLISH_MAPPED_BAM {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam, stageAs: "input/*"), path(index, stageAs: "input/*")

    output:
    tuple val(meta), path("*.bam"), path("*.{bai,csi,crai}"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}_${meta.filter_status ?: 'filter_none'}"
    """
    ln -s ${bam} ${prefix}.bam
    ln -s ${index} ${prefix}.bam.${index.getExtension()}
    """
}
