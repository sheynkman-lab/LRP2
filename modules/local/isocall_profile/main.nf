process ISOCALL_PROFILE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/pacbio/isocall:1.3.0_build1"

    input:
    tuple val(meta), path(aligned_bam)

    output:
    tuple val(meta), path("*_profile.gz"), emit: profile
    tuple val(meta), path("*_S1_PACBIO_ISOCALL_M3_ISOCALL_PROFILE_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_id = meta.sample_name ?: prefix
    def io_threads = task.cpus ?: 1
    """
    exec > >(tee ${prefix}_S1_PACBIO_ISOCALL_M3_ISOCALL_PROFILE_log.txt) 2>&1

    isocall profile \\
        --reads $aligned_bam \\
        --sample ${sample_id} \\
        --output ${prefix}_profile.gz \\
        --io-threads $io_threads \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_profile.gz
    touch ${prefix}_S1_PACBIO_ISOCALL_M3_ISOCALL_PROFILE_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
