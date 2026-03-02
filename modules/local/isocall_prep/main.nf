process ISOCALL_PREP {
    tag "$gtf.baseName"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/pacbio/isocall:0.15.0_build1' :
        'quay.io/pacbio/isocall:0.15.0_build1' }"

    input:
    path(gtf)

    output:
    path("*.isoforms.gz"), emit: isoforms
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
