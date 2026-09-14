process MEDICC2 {
    tag "$patient"
    label 'process_medium'
    maxRetries 1

    conda     "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medicc2:1.0.2--py38hcbe9525_0' :
        'biocontainers/medicc2:1.1.2--py39h0dd7abe_0' }"

    input:
    tuple val(patient), path(tsv)
    path(medicc_genes)
    path(medicc_arms)

    output:
    tuple val(patient), path("medicc2_output"),  emit: medicc2
    tuple val("${task.process}"),
      val('medicc2'),
      eval("medicc2 --version | awk 'NR == 1 { print \$NF }'"),
      emit: versions_medicc2,
      topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def plot_style = task.attempt == 1 ? 'both' : 'auto'
    """
    if [ ! -d medicc2_output ]; then
        mkdir medicc2_output
    fi

    awk '{if( \$5 < 0){sub(\$5,0)}{print}}' ${patient}.tsv > ${patient}_mod.tsv
	echo "\$?"

	medicc2 \\
    --events \\
    --chromosomes-bed $medicc_arms \\
    --regions-bed $medicc_genes \\
    --plot ${plot_style} \\
    --n-cores 4 \\
    --total-copy-numbers \\
    --input-allele-columns Copies \\
    --normal-name Diploid \\
    ${patient}_mod.tsv medicc2_output
    """
    stub:
    """
    mkdir -p medicc2_output
    """
}
