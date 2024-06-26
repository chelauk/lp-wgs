process PREP_ASCAT {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "r-ace.sif"

    input:
    tuple val(meta), path(bam), path(bai)
    val(bin)

    output:
    //tuple val(meta), path("1000kbp"), path("500kbp"), path("100kbp"),  emit: ace
    tuple val(meta), path("*.pdf"), path("*txt"), emit: qdnaseq_out
    tuple val(meta), path("*cna_segments.txt"),  path("*_bins.txt"),  emit: for_ascat
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # prefix=$prefix

    #    if [ ! -e \$prefix.bam ]; then
    #        mv $bam \$prefix.bam
    #        mv $bai \$prefix.bai
    #    fi

    QDNAseq.R ${meta.patient} ${meta.sample} $bin $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ace: 3.16
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
