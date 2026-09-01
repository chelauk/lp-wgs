process ICHORCNA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f07cec06705b4443052d3d7eaccebdbd0078366f7d074bfd4a6893980c6e2c4b/data' :
        'community.wave.seqera.io/library/r-ichorcna:0.5.1--eed4be826f05c9d4' }"

    input:
    tuple val(meta), path(wig)
    path gc_wig
    path map_wig
    path normal_wig
    path normal_background
    path centromere
    path rep_time_wig
    path exons

    output:
    tuple val(meta), path("${prefix}.RData")             , emit: rdata
    tuple val(meta), path("${prefix}.seg")               , emit: seg
    tuple val(meta), path("${prefix}.cna.seg")           , emit: cna_seg
    tuple val(meta), path("${prefix}.seg.txt")           , emit: seg_txt
    tuple val(meta), path("${prefix}.correctedDepth.txt"), emit: corrected_depth
    tuple val(meta), path("${prefix}.params.txt")        , emit: ichorcna_params
    tuple val(meta), path("${prefix}")                   , emit: output_dir
    tuple val(meta), path("${prefix}/*.pdf")             , emit: plots
    tuple val(meta), path("**/${prefix}_genomeWide.pdf") , emit: genome_plot
    path "versions.yml"                                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args       ?: ''
    prefix = task.ext.prefix       ?: "${meta.id}"
    def norm   = normal_wig        ? "normal_wig='${normal_wig}',"          : 'normal_wig=NULL,'
    def pon    = normal_background ? "normal_panel='${normal_background}'," : 'normal_panel=NULL,'
    def map    = map_wig           ? "mapWig='${map_wig}',"                 : 'mapWig=NULL,'
    def centro = centromere        ? "centromere='${centromere}',"          : ''
    def rep    = rep_time_wig      ? "repTimeWig='${rep_time_wig}',"        : 'repTimeWig=NULL,'
    def exon   = exons             ? "exons.bed='${exons}',"                : ''
    """
    #!/usr/bin/env Rscript
    library("ichorCNA")
    library("yaml")

    # Preserve the original plotting function
    original_plotCorrectionGenomeWide <- get(
        "plotCorrectionGenomeWide",
        envir = asNamespace("ichorCNA")
    )
    
    # Skip chromosome correction plots when median read count is zero
    safe_plotCorrectionGenomeWide <- function(
        correctOutput,
        chr = NULL,
        seqinfo = NULL,
        ...
    ) {
        x <- correctOutput
    
        if (!is.null(chr)) {
            x <- x[
                as.character(seqnames(x)) ==
                    as.character(chr)
            ]
        }
    
        med <- median(x\$reads, na.rm = TRUE)
    
        if (
            length(x) == 0 ||
            !is.finite(med) ||
            med <= 0
        ) {
            warning(
                "Skipping correction plot for ",
                if (is.null(chr)) "genome-wide data" else chr,
                ": median read count is zero or unavailable"
            )
    
            return(invisible(NULL))
        }
    
        original_plotCorrectionGenomeWide(
            correctOutput = correctOutput,
            chr = chr,
            seqinfo = seqinfo,
            ...
        )
    }
    
    # Override only for this R session
    assignInNamespace(
        "plotCorrectionGenomeWide",
        safe_plotCorrectionGenomeWide,
        ns = "ichorCNA"
    )


    run_ichorCNA(
        tumor_wig='${wig}',
        id='${prefix}',
        cores=${task.cpus},
        gcWig='${gc_wig}',
        $norm
        $pon
        $map
        $centro
        $rep
        $exon
        $args
        outDir="."
    )


    ### Make Versions YAML for NF-Core ###
    versions = list()
    versions["r"]        <- paste(R.Version()\$major, R.Version()\$minor, sep=".")
    versions["ichorCNA"] <- paste(packageVersion("ichorCNA"), sep=".")

    yaml_str <- as.yaml(
        list(
            "${task.process}" = versions
        )
    )
    writeLines(yaml_str, file("versions.yml"))
    """

    stub:
    prefix = task.ext.prefix   ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library("ichorCNA")
    library("yaml")

    file.create('${prefix}.RData')
    file.create('${prefix}.seg')
    file.create('${prefix}.cna.seg')
    file.create('${prefix}.seg.txt')
    file.create('${prefix}.correctedDepth.txt')
    file.create('${prefix}.params.txt')
    dir.create('${prefix}')
    file.create('${prefix}/${prefix}_genomeWide.pdf')

    ### Make Versions YAML for NF-Core ###
    versions = list()
    versions["r"]        <- paste(R.Version()\$major, R.Version()\$minor, sep=".")
    versions["ichorCNA"] <- paste(packageVersion("ichorCNA"), sep=".")

    yaml_str <- as.yaml(
        list(
            "${task.process}" = versions
        )
    )
    writeLines(yaml_str, file("versions.yml"))
    """

}
