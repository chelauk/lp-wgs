process PREP_MEDICC2 {
    tag "$patient"
    label 'process_medium'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "docker.io/chelauk/ace-qdnaseq:1.0.1"

    input:
    tuple val(patient), val(samples), val(ids), val(ploidy),  path(ace_out)
    path(bin_dir)
    val(bin)

    output:
    tuple val(patient), path("*tsv"), emit: for_medicc
    tuple val(patient), path("*txt"), emit: for_report
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    prep_medicc.R $patient "${ploidy.join(' ')}" $bin_dir $bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medicc2: 3.16
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    """
    echo "prep_medicc.R $patient $ploidy $bin_dir $bin"

    touch ${patient}.tsv
    touch ${patient}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medicc2: stub version
    END_VERSIONS
    """
}
