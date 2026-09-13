process RUN_BAYES {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/chelauk/bcp-qdnaseq:1.0.0"

    input:
    tuple val(meta), path(bam), path(bai)
    val (bin)
    val(qdnaseq_genome)
    path(bin_dir)

    output:
    tuple val(meta), path("${meta.patient}_${meta.sample}_bcp"), emit: bayes_cn
    tuple val("${task.process}"),
          val('r-base'),
          eval("Rscript --vanilla -e 'cat(as.character(getRversion()))'"),
          emit: versions_r,
          topic: versions
    
    tuple val("${task.process}"),
          val('bcp'),
          eval("Rscript --vanilla -e 'cat(as.character(packageVersion(\"bcp\")))'"),
          emit: versions_bcp,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome = qdnaseq_genome ?: 'hg38'
    """
    bayescnasketch.R ${meta.patient} ${meta.sample} $bin ${bin_dir} $bam ${genome}
    """

    stub:
    def genome = qdnaseq_genome ?: 'hg38'
    """
    echo  "bayescnasketch.R ${meta.patient} ${meta.sample} $bin ${bin_dir} ${genome}"
    mkdir -p ${meta.patient}_${meta.sample}_bcp
    touch ${meta.patient}_${meta.sample}_bcp/${meta.patient}_${meta.sample}_bcp_segments.csv
    touch ${meta.patient}_${meta.sample}_bcp/${meta.patient}_${meta.sample}_bcp_wgs_profile_selected.pdf
    """
}
