process ISOCALL_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/pacbio/isocall:0.15.0_build1' :
        'quay.io/pacbio/isocall:0.15.0_build1' }"

    input:
    tuple val(meta), path(profiles)

    output:
    tuple val(meta), path("*_merged_profiles.gz"), emit: merged_profile
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def profile_list = profiles.collect{ it }.join(' ')
    """
    isocall merge \\
        --profiles $profile_list \\
        --output ${prefix}_merged_profiles.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged_profiles.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
