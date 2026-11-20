process ANEUFINDER {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/bioconductor-aneufinder:1.34.0--r44he5774e6_0'

    input:
    tuple val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*')
    val(bin_size_kb)
    val(paired_end)

    output:
    tuple val(meta), path('aneufinder'), emit: results
    tuple val("${task.process}"), val('aneufinder'),
        eval('Rscript --vanilla -e \'cat(as.character(packageVersion("AneuFinder")))\''),
        emit: versions_aneufinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bin_size_bp = bin_size_kb.toInteger() * 1000
    def paired_end_r = paired_end ? 'TRUE' : 'FALSE'
    """
    Rscript --vanilla - <<'RSCRIPT'
    set.seed(1)
    AneuFinder::Aneufinder(
        inputfolder = "input",
        outputfolder = "aneufinder",
        binsizes = ${bin_size_bp},
        pairedEndReads = ${paired_end_r},
        numCPU = ${task.cpus},
        reuse.existing.files = FALSE,
        hotspot.pval = NULL,
        cluster.plots = FALSE
    )
    RSCRIPT
    """

    stub:
    """
    mkdir -p aneufinder
    """
}
