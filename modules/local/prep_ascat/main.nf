process RUN_QDNASEQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "docker.io/chelauk/ace-qdnaseq:1.0.1"

    input:
    tuple val(meta), path(bam), path(bai)
    val(bin)
    val(qdnaseq_genome)
    val(qdnaseq_package)

    output:
    tuple val(meta), path("*.pdf"), path("*.txt"), emit: qdnaseq_out
    tuple val(meta), path("*cna_segments.txt"),  path("*bins.txt"),  emit: for_ascat
    tuple val(meta), path("*.rds"), emit: for_ace
    tuple val("${task.process}"),
          val('r-base'),
          eval("Rscript --vanilla -e 'cat(as.character(getRversion()))'"),
          emit: versions_r,
          topic: versions
    
    tuple val("${task.process}"),
          val('r-qdnaseq'),
          eval("Rscript --vanilla -e 'cat(as.character(packageVersion(\"QDNAseq\")))'"),
          emit: versions_qdnaseq,
          topic: versions
    
    tuple val("${task.process}"),
          val('r-cghcall'),
          eval("Rscript --vanilla -e 'cat(as.character(packageVersion(\"CGHcall\")))'"),
          emit: versions_cghcall,
          topic: versions
    
    tuple val("${task.process}"),
          val("r-${(qdnaseq_package ?: 'QDNAseq.hg38').toLowerCase()}"),
          eval(
              "Rscript --vanilla -e 'cat(as.character(packageVersion(\"${qdnaseq_package ?: 'QDNAseq.hg38'}\")))'"
          ),
          emit: versions_qdnaseq_reference,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome = qdnaseq_genome ?: 'hg38'
    def qdnaseqPackage = qdnaseq_package ?: 'QDNAseq.hg38'
    """
    QDNAseq.R ${meta.patient} ${meta.sample} $bin $bam ${genome} ${qdnaseqPackage}
    """

    stub:
    def genome = qdnaseq_genome ?: 'hg38'
    def qdnaseqPackage = qdnaseq_package ?: 'QDNAseq.hg38'
    """
    echo  "QDNAseq.R ${meta.patient} ${meta.sample} $bin $bam ${genome} ${qdnaseqPackage}"
    touch "${meta.id}.cna_segments.txt"
    touch "${meta.id}.bins.txt"
    touch "${meta.id}_${bin}kbp.rds"
    touch "${meta.id}.called_segments.pdf"
    """
}
