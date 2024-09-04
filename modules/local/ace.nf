process ACE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "r-ace.sif"

    input:
    tuple val(meta), path(bam), path(bai)
    val(filter_status)

    output:
    //tuple val(meta), path("1000kbp"), path("500kbp"), path("100kbp"),  emit: ace
    tuple val(meta), path("${meta.sample}_${filter_status}"),  emit: ace
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prefix=$prefix

    if [ ! -e \$prefix.bam ]; then
        mv $bam \$prefix.bam
    fi

    if [ ! -d "${meta.sample}_${filter_status}" ]; then
    mkdir "${meta.sample}_${filter_status}"
    fi
    ace.R ${meta.sample}_${filter_status}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ace: 3.16
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${meta.sample}_${filter_status}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ace: stub version
    END_VERSIONS
    """
}
