process ISOCALL_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/pacbio/isocall:0.15.0_build1"

    input:
    tuple val(meta), path(merged_profile)
    path(known_isoforms)
    path(reference_fasta)
    path(config_toml)

    output:
    tuple val(meta), path("*.isocall.isoforms.gtf.gz"), emit: gtf
    tuple val(meta), path("*.isocall.count_matrix.txt"), emit: count_matrix
    tuple val(meta), path("*_S1_PACBIO_ISOCALL_M5_ISOCALL_CALL_log.txt"), emit: log
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
    exec > >(tee ${prefix}_S1_PACBIO_ISOCALL_M5_ISOCALL_CALL_log.txt) 2>&1

    isocall call \\
        --threads $threads \\
        --merged-profile $merged_profile \\
        --known-isoforms $known_isoforms \\
        --reference $reference_fasta \\
        --output-prefix ${prefix}.isocall \\
        --config $config_toml \\
        --min-reads-per-isoform $min_read_support \\
        --max-bundles-per-gene $max_bundles_per_gene \\
        $args
    
    mv ${prefix}.isocall.count_matrix.txt ${prefix}.isocall.count_matrix.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.isocall.isoforms.gtf.gz
    touch ${prefix}.isocall.count_matrix.txt
    touch ${prefix}_S1_PACBIO_ISOCALL_M5_ISOCALL_CALL_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
