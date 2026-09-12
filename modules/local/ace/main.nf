process ACE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconductor-ace:1.16.0--r42hdfd78af_0"
    container "docker.io/chelauk/ace-qdnaseq:1.0.1"

    input:
    tuple val(meta), path(qdnaseq_rds)
    val(filter_status)
    val(qdnaseq_genome)
    val(ploidy)
    val(bin)

    output:
    tuple val(meta), path("${meta.sample}_${filter_status}"),  emit: ace
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome = qdnaseq_genome ?: 'hg38'
    """
    if [ ! -d "${meta.sample}_${filter_status}" ]; then
    mkdir "${meta.sample}_${filter_status}"
    fi
    ace.R ${meta.sample}_${filter_status} ${genome} "${ploidy}" ${bin} .
    if [ ! -e "${meta.sample}_${filter_status}/${bin}kbp.rds" ]; then
        cp ${qdnaseq_rds} "${meta.sample}_${filter_status}/${bin}kbp.rds"
    fi
    
    Rscript -e '
     versions <- list(
         "r-base" = as.character(getRversion()),
         "r-ACE"  = as.character(packageVersion("ACE"))
     )
     yaml::write_yaml(
         setNames(list(versions), "${task.process}"),
         "versions.yml"
     )
    ' 
    """
    stub:
    def genome = qdnaseq_genome ?: 'hg38'
    """
    mkdir ${meta.sample}_${filter_status}
    echo "ace.R ${meta.sample}_${filter_status} ${genome} \"${ploidy}\" ${bin} ."
    Rscript -e '
     versions <- list(
         "r-base" = as.character(getRversion()),
         "r-ACE"  = as.character(packageVersion("ACE"))
     )
     yaml::write_yaml(
         setNames(list(versions), "${task.process}"),
         "versions.yml"
     )
    ' 
    """
}
