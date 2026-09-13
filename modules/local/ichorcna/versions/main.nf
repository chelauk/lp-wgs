process ICHORCNA_VERSIONS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f07cec06705b4443052d3d7eaccebdbd0078366f7d074bfd4a6893980c6e2c4b/data' :
        'community.wave.seqera.io/library/r-ichorcna:0.5.1--eed4be826f05c9d4' }"

    output:
    tuple val("${task.process}"),
          val('r-base'),
          eval("Rscript --vanilla -e 'cat(as.character(getRversion()))'"),
          topic: versions,
          emit: versions_r

    tuple val("${task.process}"),
          val('r-ichorCNA'),
          eval("Rscript --vanilla -e 'cat(as.character(packageVersion(\"ichorCNA\")))'"),
          topic: versions,
          emit: versions_ichorcna

    script:
    """
    true
    """
}
