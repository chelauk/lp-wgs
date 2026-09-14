process SAMBAMBA_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda     "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:0.8.1--hadffe2f_1' :
        'quay.io/biocontainers/sambamba:0.8.1--hadffe2f_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    tuple val("${task.process}"), val('sambamba'), eval("sambamba --version 2>&1 | grep -oPm1 'sambamba \\\\K[0-9.]+'"), topic: versions, emit: versions_sambamba

    script:
    """
    sambamba merge \\
    --nthreads=4 /dev/stdout ${bam.join(' ')} | \\
    sambamba sort --tmpdir . -o ${meta.id}.bam /dev/stdin
    """

    stub:
    """
    echo -e "sambamba merge \\
    --nthreads=4 /dev/stdout ${bam.join(' ')} | \\
    sambamba sort --tmpdir . -o ${meta.id}.bam /dev/stdin"
    touch ${meta.id}.bam
    """
}
