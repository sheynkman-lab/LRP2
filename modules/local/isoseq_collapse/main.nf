process ISOSEQ_COLLAPSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:4.0.0--h9ee0642_0' :
        'biocontainers/isoseq:4.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(aligned_bam), path(flnc_bam), path(aligned_bam_pbi), path(flnc_bam_pbi)

    output:
    tuple val(meta), path("*.collapsed.gff"), emit: gff
    tuple val(meta), path("*.collapsed.flnc_count.txt"), emit: flnc_count
    tuple val(meta), path("*.collapsed.group.txt"), emit: group_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def max_fuzzy_junction = task.ext.max_fuzzy_junction ?: params.max_fuzzy_junction
    def max_5p_diff = task.ext.max_5p_diff ?: params.max_5p_diff
    def max_3p_diff = task.ext.max_3p_diff ?: params.max_3p_diff
    """
    isoseq collapse \\
        --max-fuzzy-junction $max_fuzzy_junction \\
        --max-5p-diff $max_5p_diff \\
        --max-3p-diff $max_3p_diff \\
        $aligned_bam \\
        $flnc_bam \\
        ${prefix}.collapsed.gff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$( isoseq collapse --version | head -n 1 | sed 's/isoseq collapse //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.collapsed.gff
    touch ${prefix}.collapsed.flnc_count.txt
    touch ${prefix}.collapsed.group.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$( isoseq collapse --version | head -n 1 | sed 's/isoseq collapse //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}

