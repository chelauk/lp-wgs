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

    tuple val("${task.process}"),
      val('forecast-ascat'),
      eval('forecast_ascat.R --version'),
      emit: versions_forecast_ascat,
      topic: versions

    tuple val("${task.process}"),
      val('r-base'),
      eval("Rscript --vanilla -e 'cat(as.character(getRversion()))'"),
      emit: versions_r,
      topic: versions

    tuple val("${task.process}"),
      val('r-copynumber'),
      eval("Rscript --vanilla -e 'cat(as.character(packageVersion(\"copynumber\")))'"),
      emit: versions_copynumber,
      topic: versions

    tuple val("${task.process}"),
      val('r-ggplot2'),
      eval("Rscript --vanilla -e 'cat(as.character(packageVersion(\"ggplot2\")))'"),
      emit: versions_ggplot2,
      topic: versions

    tuple val("${task.process}"),
      val('r-cowplot'),
      eval("Rscript --vanilla -e 'cat(as.character(packageVersion(\"cowplot\")))'"),
      emit: versions_cowplot,
      topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome = qdnaseq_genome ?: 'hg38'

    def ploidyCommands = ploidies.toString().split(',')  
        .collect { p -> p.trim() }  
        .findAll { p -> p }  
        .collect { ploidy ->  
        def ploidyDir = "ascat_ploidy_${ploidy.replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
    """
    mkdir -p "${ploidyDir}"
    cp ${cna_segments} ${cna_bins} "${ploidyDir}/"
    (
        cd "${ploidyDir}"
        forecast_ascat.R \
            ${meta.id} \
            ${ploidy} \
            \$(dirname "\$(command -v forecast_ascat.R)") \
            ../${chr_arm_boundaries} \
            ${genome} \
            ${ascat_pcf_gamma}
    )
    """.stripIndent().trim()   
     }.join('\n')
    """
    ${ploidyCommands}

    """

    stub:
    def genome = qdnaseq_genome ?: 'hg38'
    def ploidyCommands = ploidies.toString().split(',')  
        .collect { p -> p.trim() }  
        .findAll { p -> p }  
        .collect { ploidy ->  
        def ploidyDir = "ascat_ploidy_${ploidy.replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
        """
        mkdir -p "${ploidyDir}"
        echo "forecast_ascat.R ${meta.id} ${ploidy} ../${chr_arm_boundaries} ${genome} ${ascat_pcf_gamma}" > "${ploidyDir}/${meta.id}_ploidy_${ploidy}.command.txt"
        touch "${ploidyDir}/${meta.id}_selected_ascat_lp_plot.pdf"
        touch "${ploidyDir}/${meta.id}_selected_cna_ploidy_search_calls.txt"
        touch "${ploidyDir}/${meta.id}_selected_ascat_lp_metrics.txt"
        """.stripIndent().trim()
    }.join('\n')
    """
    ${ploidyCommands}

    """
}
