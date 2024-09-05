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
    tuple val(meta), path("*cna_segments.txt"),  path("*bins.txt"),  emit: for_ascat
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    QDNAseq.R ${meta.patient} ${meta.sample} $bin $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prep_ascat: 1
    END_VERSIONS
    """
    
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${meta.id}.cna_segments.txt"
    touch "${meta.id}.bins.txt"
    touch "${meta.id}.called_segments.pdf"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prep_ascat: stub version
    END_VERSIONS
    """
}
