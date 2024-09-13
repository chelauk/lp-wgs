process RUN_ASCAT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "forecast_lp_ascat.sif"
    
    input:
    tuple val(meta), path(cna_segments), path(cna_bins)
    val(ploidy)
	path(chr_arm_boundaries)

    output:
    tuple val(meta), path("*txt"), path("*pdf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    forecast_ascat.R ${meta.id} ${ploidy} ${projectDir}/bin/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: 0.01
    END_VERSIONS
    """
    
    stub:
    def args = task.ext.args ?: ''
    """
    echo  "forecast_ascat.R ${meta.id} ${ploidy} ${chr_arm_boundaries}"
    touch ${meta.id}.pdf
    touch ${meta.id}.txt    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: stub version
    END_VERSIONS
    """
}
