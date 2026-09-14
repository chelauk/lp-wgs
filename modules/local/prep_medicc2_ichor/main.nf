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
    tuple val("${task.process}"),
      val('prep-medicc-ichor'),
      eval("python3 ${bin_dir}/prep_medicc_ichor.py --version"),
      emit: versions_prep_medicc_ichor,
      topic: versions

    tuple val("${task.process}"),
      val('python'),
      eval("python3 --version | sed 's/Python //'"),
      emit: versions_python,
      topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${bin_dir}/prep_medicc_ichor.py \\
        --patient ${patient} \\
        --out ${patient}.tsv \\
        --report medicc2_ichor_prep.txt \\
        ${segs}
    """

    stub:
    """
    touch ${patient}.tsv
    touch medicc2_ichor_prep.txt
    """
}
