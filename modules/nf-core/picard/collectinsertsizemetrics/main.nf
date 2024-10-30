process PICARD_COLLECTINSERTSIZEMETRICS {
    tag "$meta.id"
    label "process_long"

    conda "bioconda::picard=2.26.10" 
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val filter_status

    output:
    tuple val(meta), path("*.txt"),  emit: size_metrics
    tuple val(meta), path("*.pdf"), emit: size_plots
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bams
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3
    if (!task.memory) {
        log.info "[Picard CollectInsertSizeMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
    } else {
        avail_mem = task.memory.giga
    }
    """
    
    bam_link=\$( readlink $bam )
    ln -s \$bam_link ${prefix}_${filter_status}_${args}.bam 
    bai_link=\$( readlink $bai )
    ln -s \$bai_link ${prefix}_${filter_status}_${args}.bam.bai 
    
    picard \\
    CollectInsertSizeMetrics \\
    I=$bam \\
    O=${prefix}_${filter_status}_${args}.insert_sizes.txt \\
    H=${prefix}_${filter_status}_${args}.insert_sizes.pdf \\
    HISTOGRAM_WIDTH=${args} \\
    M=0.5 \\
    VALIDATION_STRINGENCY=SILENT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectInsertSizeMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ln -s $bam ${prefix}_${filter_status}_processed.bam
    ln -s $bai ${prefix}_${filter_status}_processed.bai
    echo ${args}
    touch ${prefix}_${filter_status}_${args}.collectinsertsizemetrics.txt
    touch ${prefix}_${filter_status}_${args}.collectinsertsizemetrics.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 2.26.10 
    END_VERSIONS
    """
}
