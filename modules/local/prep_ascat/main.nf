process PREP_ASCAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "ace-qdnaseq_latest.sif"

    input:
    tuple val(meta), path(bam), path(bai)
    val(bin)
    val(qdnaseq_genome)
    val(qdnaseq_package)

    output:
    tuple val(meta), path("*.pdf"), path("*txt"), emit: qdnaseq_out
    tuple val(meta), path("*cna_segments.txt"),  path("*bins.txt"),  emit: for_ascat
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    def qdnaseqPackage = qdnaseq_package ?: 'QDNAseq.hg38'
    """
    QDNAseq.R ${meta.patient} ${meta.sample} $bin $bam ${genome} ${qdnaseqPackage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prep_ascat: 1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def genome = qdnaseq_genome ?: 'hg38'
    def qdnaseqPackage = qdnaseq_package ?: 'QDNAseq.hg38'
    """
    echo  "QDNAseq.R ${meta.patient} ${meta.sample} $bin $bam ${genome} ${qdnaseqPackage}"
    touch "${meta.id}.cna_segments.txt"
    touch "${meta.id}.bins.txt"
    touch "${meta.id}.called_segments.pdf"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prep_ascat: stub version
    END_VERSIONS
    """
}
