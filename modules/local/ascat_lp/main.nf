process RUN_ASCAT {
    tag "$patient"
    label 'process_medium'

    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "ascatlp.0.1_latest.sif"
    
    input:
    tuple val(patient), val(samples), val(ids), path(cna_segments), path(cna_bins)
    path(chr_arm_boundaries)

    output:
    tuple val(patient), val(samples), path("*txt"), path("*pdf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    forecast_ascat.R $patient "${samples.join(' ')}" "${ids.join(' ')}" ${chr_arm_boundaries}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: 0.01
    END_VERSIONS
    """
    
    stub:
    def args = task.ext.args ?: ''
    """
    touch ${patient}.pdf
    touch ${patient}.txt    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: stub version
    END_VERSIONS
    """
}
