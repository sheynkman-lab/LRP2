process ISOCALL_PREP {
    tag "$gtf.baseName"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/pacbio/isocall:0.15.0_build1"

    input:
    path(gtf)

    output:
    path("*.isoforms.gz"), emit: isoforms
    path("*_S1_PACBIO_ISOCALL_M1_ISOCALL_PREP_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    exec > >(tee ${prefix}_S1_PACBIO_ISOCALL_M1_ISOCALL_PREP_log.txt) 2>&1

    isocall prep-isoforms \\
        --gtf $gtf \\
        --output ${prefix}.isoforms.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}.isoforms.gz
    touch ${prefix}_S1_PACBIO_ISOCALL_M1_ISOCALL_PREP_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
