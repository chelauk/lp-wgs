process RUN_ASCAT {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/chelauk/ascat_lp:1.0.1"

    input:
    tuple val(meta), path(cna_segments), path(cna_bins)
    val(ploidy)
    path(chr_arm_boundaries)
    val(qdnaseq_genome)
    val(ascat_pcf_gamma)

    output:
    tuple val(meta), path("*txt"), path("*pdf")
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    """
    forecast_ascat.R ${meta.id} ${ploidy} ${projectDir}/bin/ ${chr_arm_boundaries} ${genome} ${ascat_pcf_gamma}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: 0.01
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    """
    echo  "forecast_ascat.R ${meta.id} ${ploidy} ${projectDir}/bin/ ${chr_arm_boundaries} ${genome} ${ascat_pcf_gamma}"
    touch ${meta.id}.pdf
    touch ${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: stub version
    END_VERSIONS
    """
}
