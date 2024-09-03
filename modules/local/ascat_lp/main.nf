process RUN_ASCAT {
    tag "$patient"
    label 'process_medium'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "ascatlp.0.1_latest.sif"

    input:
    tuple val(patient), val(samples), val(ids), path(cna_segments)
    val(ploidy)
    val(purity)


    output:
    tuple val(patient), val(samples), path("*rds"), path("*txt"), path("*pdf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ascat_lp.R $patient "${samples.join(' ')}" "${ids.join(' ')}" $ploidy $purity  ${workflow.projectDir}/bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ascat_lp: 0.01
    END_VERSIONS
    """
    
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $filter_status

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ace: stub version
    END_VERSIONS
    """
}
