process MEDICC2 {
    tag "$patient"
    label 'process_medium'
    maxRetries 1

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconda::medicc2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medicc2:1.0.2--py38hcbe9525_0' :
        'biocontainers/medicc2:1.0.2--py38hcbe9525_0' }"

    input:
    tuple val(patient), path(tsv)

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
    sed -i 's/-[0-9]*/0/' ${patient}.tsv
    medicc2 \\
    --plot ${plot_style} \\
    --n-cores 4 \\
    --input-allele-columns Copies \\
    ${patient}.tsv medicc2_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medicc2: \$(medicc2 --version | sed -n 1p | cut -d ' ' -f 1,2)
    END_VERSIONS
    """
    stub:
    def plot_style = task.attempt == 1 ? 'heatmap' : 'auto'
    def args = task.ext.args ?: ''
    """
    if [ ! -d medicc2_output ]; then
        mkdir medicc2_output
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medicc2: stub version
    END_VERSIONS
    """
}
