process PICARD_COLLECTALIGNMENTSUMMARYMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::picard=3.0.0 r::r-base"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path dict

    output:
    tuple val(meta), path("*metrics"), emit: metrics
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectAlignmentSummaryMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectAlignmentSummaryMetrics \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.CollectAlignmentSummaryMetrics.metrics \\
        --REFERENCE_SEQUENCE ${fasta} 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectAlignmentSummaryMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectAlignmentSummaryMetrics.metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectAlignmentSummaryMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
