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
    tuple val(meta), path("ichor_*")   , emit: ichor_out
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? "${task.ext.args} ${normal_wig}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    runIchorCNA.R \\
        $args \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir ./ichor_"${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ? "${task.ext.args} ${normal_wig}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo -e "runIchorCNA.R \\
        $args \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir . 

    cp */*genomeWide.pdf ."
    touch ${prefix}.cna.seg
    touch ${prefix}.params.txt
    touch ${prefix}genomeWide.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
}
