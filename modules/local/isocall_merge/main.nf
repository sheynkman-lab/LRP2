process ISOCALL_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/pacbio/isocall:1.3.0_build1"

    input:
    tuple val(meta), path(profiles)

    output:
    tuple val(meta), path("*.isocall.merged_profiles.gz"), emit: merged_profile
    tuple val(meta), path("*_S1_PACBIO_ISOCALL_M4_ISOCALL_MERGE_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def profile_list = profiles.collect{ it }.join(' ')
    """
    exec > >(tee ${prefix}_S1_PACBIO_ISOCALL_M4_ISOCALL_MERGE_log.txt) 2>&1

    isocall merge \\
        --profiles $profile_list \\
        --output ${prefix}.isocall.merged_profiles.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.isocall.merged_profiles.gz
    touch ${prefix}_S1_PACBIO_ISOCALL_M4_ISOCALL_MERGE_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
