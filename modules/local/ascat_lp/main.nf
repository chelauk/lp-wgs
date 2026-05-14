process RUN_ASCAT {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/chelauk/ascat_lp:1.0.1"

    input:
    tuple val(meta), path(cna_segments), path(cna_bins)
    val(ploidies)
    path(chr_arm_boundaries)
    val(qdnaseq_genome)
    val(ascat_pcf_gamma)

    output:
    tuple val(meta), path("ascat_ploidy_*"), emit: ascat
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    def ploidyCommands = ploidies.toString().split(',').collect { it.trim() }.findAll { it }.collect { ploidy ->
        def ploidyDir = "ascat_ploidy_${ploidy.replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
        """
        mkdir -p "${ploidyDir}"
        cp ${cna_segments} ${cna_bins} "${ploidyDir}/"
        (
            cd "${ploidyDir}"
            forecast_ascat.R ${meta.id} ${ploidy} ${projectDir}/bin/ ../${chr_arm_boundaries} ${genome} ${ascat_pcf_gamma}
        )
        """.stripIndent().trim()
    }.join('\n')
    """
    ${ploidyCommands}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: 0.01
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    def ploidyCommands = ploidies.toString().split(',').collect { it.trim() }.findAll { it }.collect { ploidy ->
        def ploidyDir = "ascat_ploidy_${ploidy.replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
        """
        mkdir -p "${ploidyDir}"
        echo "forecast_ascat.R ${meta.id} ${ploidy} ${projectDir}/bin/ ../${chr_arm_boundaries} ${genome} ${ascat_pcf_gamma}" > "${ploidyDir}/${meta.id}_ploidy_${ploidy}.command.txt"
        touch "${ploidyDir}/${meta.id}_selected_ascat_lp_plot.pdf"
        touch "${ploidyDir}/${meta.id}_selected_cna_ploidy_search_calls.txt"
        touch "${ploidyDir}/${meta.id}_selected_ascat_lp_metrics.txt"
        """.stripIndent().trim()
    }.join('\n')
    """
    ${ploidyCommands}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: stub version
    END_VERSIONS
    """
}
