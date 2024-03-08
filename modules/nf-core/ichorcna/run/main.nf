process ICHORCNA_RUN {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::r-ichorcna=0.3.2"
    container "ichorcna_1_1.sif"

    input:
    tuple val(meta), path(wig)
    path normal_wig
    path(gc_wig)
    path map_wig
    path panel_of_normals
    path centromere
    val filter_status

    output:
    tuple val(meta), path("filter*"),  emit: ichor_out
    path "versions.yml"             ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? "${task.ext.args} ${normal_wig}" : ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ ! -d $filter_status ]; then
    mkdir -p $filter_status
    fi
    runIchorCNA.R \\
        $args \\
        $args2 \\
        $args3 \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        --genomeBuild "hg38" \\
        ${pon} \\
        ${centro} \\
        --outDir $filter_status

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ? "${task.ext.args} ${normal_wig}" : ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ ! -d $filter_status ]; then
    mkdir -p $filter_status
    fi
    echo -e 'runIchorCNA.R \\
        $args \\
        $args2 \\
        $args3 \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        --outDir ${filter_status}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
}
