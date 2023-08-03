process ICHORCNA_RUN {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::r-ichorcna=0.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.3.2--pl5321r42hdfd78af_2' :
        'quay.io/biocontainers/r-ichorcna:0.3.2--pl5321r42hdfd78af_2' }"

    input:
    tuple val(meta), path(wig)
    path normal_wig
    path(gc_wig)
    path map_wig
    path panel_of_normals
    path centromere

    output:
    tuple val(meta), path("solutions"),  emit: ichor_out
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? "${task.ext.args} ${normal_wig}" : ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ ! -d solutions ]; then
    mkdir -p solutions
    fi
    runIchorCNA.R \\
        $args \\
        $args2 \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir solutions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ? "${task.ext.args} ${normal_wig}" : ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo -e "runIchorCNA.R \\
        $args \\
        $args2 \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir solutions"
    if [ ! -d solutions ] 
        then
            mkdir -p solutions
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
}
