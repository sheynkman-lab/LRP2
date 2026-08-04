process SANITIZE_PAR_IDS {
    tag "${gtf}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*_sanitized.gtf"), emit: sanitized_gtf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output_name = "${gtf.baseName}_sanitized.gtf"
    """
    # Replace _PAR_Y with -PAR-Y to fix issue with PAR gene IDs (affects GENCODE references v25 - v43)
    sed 's/_PAR_Y/-PAR-Y/g' ${gtf} > ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | head -n1 | sed 's/.*) //')
    END_VERSIONS
    """

    stub:
    def output_name = "${gtf.baseName}_sanitized.gtf"
    """
    touch ${output_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | head -n1 | sed 's/.*) //')
    END_VERSIONS
    """
}
