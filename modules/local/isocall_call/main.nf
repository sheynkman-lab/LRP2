process ISOCALL_CALL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/pacbio/isocall:0.15.0_build1' :
        'quay.io/pacbio/isocall:0.15.0_build1' }"

    input:
    tuple val(meta), path(merged_profile)
    path(known_isoforms)
    path(reference_fasta)
    path(config_toml)

    output:
    tuple val(meta), path("*.isoforms.gtf.gz"), emit: gtf
    tuple val(meta), path("*.count_matrix.txt"), emit: count_matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 1
    def min_read_support = task.ext.min_read_support ?: params.min_read_support
    def max_bundles_per_gene = task.ext.max_bundles_per_gene ?: params.max_bundles_per_gene
    """
    isocall call \\
        --threads $threads \\
        --merged-profile $merged_profile \\
        --known-isoforms $known_isoforms \\
        --reference $reference_fasta \\
        --output-prefix ${prefix} \\
        --config $config_toml \\
        --min-reads-per-isoform $min_read_support \\
        --max-bundles-per-gene $max_bundles_per_gene \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.isoforms.gtf.gz
    touch ${prefix}.count_matrix.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
