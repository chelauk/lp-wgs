process RUN_BAYES {
    tag "${meta.id}"
    label 'process_low'

    container "bcp-qdnaseq_latest.sif"

    input:
    tuple val(meta), path(bam), path(bai)
    val (bin)
    val(qdnaseq_genome)

    output:
    tuple val(meta), path("${meta.patient}_${meta.sample}_bcp"), emit: bayes_cn
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    """
    bayescnasketch.R ${meta.patient} ${meta.sample} $bin ${projectDir}/bin/ $bam ${genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcp: 0.01
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    """
    echo  "bayescnasketch.R ${meta.patient} ${meta.sample} $bam ${projectDir}/bin/ ${genome}"
    mkdir -p ${meta.patient}_${meta.sample}_bcp
    touch ${meta.patient}_${meta.sample}_bcp/${meta.patient}_${meta.sample}_bcp_segments.csv
    touch ${meta.patient}_${meta.sample}_bcp/${meta.patient}_${meta.sample}_bcp_wgs_profile_selected.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcp: stub version
    END_VERSIONS
    """
}
