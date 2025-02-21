process MEDICC2 {
    tag "$patient"
    label 'process_medium'
    maxRetries 1

    conda '/home/chela.james/miniconda3/envs/medicc2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medicc2:1.0.2--py38hcbe9525_0' :
        'biocontainers/medicc2:1.1.2--py39h0dd7abe_0' }"

    input:
    tuple val(patient), path(tsv)
    path(medicc_genes)
    path(medicc_arms)

    output:
    tuple val(patient), path("medicc2_output"),  emit: medicc2
    path "versions.yml"                    ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def plot_style = task.attempt == 1 ? 'both' : 'auto'
    def args = task.ext.args ?: ''
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medicc2: \$(medicc2 --version | sed -n 1p | cut -d ' ' -f 1,2)
    END_VERSIONS
    """
    stub:
    def plot_style = task.attempt == 1 ? 'heatmap' : 'auto'
    def args = task.ext.args ?: ''
    """
    mkdir -p medicc2_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    medicc2: stub version
    END_VERSIONS
    """
}
