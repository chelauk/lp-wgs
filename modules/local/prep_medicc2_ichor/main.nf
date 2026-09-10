process PREP_MEDICC2_ICHOR {
    tag "$patient"
    label 'process_low'

    container "python:3.11-bookworm"

    input:
    tuple val(patient), val(samples), val(ids), path(segs)
    path(bin_dir)

    output:
    tuple val(patient), path("${patient}.tsv"), emit: for_medicc
    tuple val(patient), path("medicc2_ichor_prep.txt"), emit: for_report
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${bin_dir}prep_medicc_ichor.py \\
        --patient ${patient} \\
        --out ${patient}.tsv \\
        --report medicc2_ichor_prep.txt \\
        ${args} \\
        ${segs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${patient}.tsv
    touch medicc2_ichor_prep.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub version
    END_VERSIONS
    """
}
