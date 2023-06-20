process ACE {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-ace:1.16.0--r42hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-ace:1.16.0--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}"), emit: ace
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prefix=$prefix
    if [ -d \$prefix ]
    then
        echo "Directory /path/to/dir exists."
    else
        echo "Error: Directory /path/to/dir does not exist."
        mkdir \$prefix
    fi

    ace.R \$prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ace: 3.16
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prefix=$prefix
    mkdir \$prefix
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ace: stub version
    END_VERSIONS
    """
}
