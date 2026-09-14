process ISOCALL_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jtllab/isocall:1.3.0-nextflow-fix' :
        'jtllab/isocall:1.3.0-nextflow-fix' }"

    input:
    tuple val(meta), path(merged_profile)
    path(known_isoforms)
    path(reference_fasta)
    path(config_toml)

    output:
    tuple val(meta), path("*.isocall.isoforms.gtf.gz"), emit: gtf
    tuple val(meta), path("*.isocall.count_matrix.csv"), emit: count_matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 1
    """
    isocall call \\
        --threads $threads \\
        --merged-profile $merged_profile \\
        --known-isoforms $known_isoforms \\
        --reference $reference_fasta \\
        --output-prefix ${prefix}.isocall \\
        --config $config_toml \\
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
    touch ${prefix}.isocall.count_matrix.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isocall: \$( isocall --version 2>&1 | sed 's/isocall //g' )
    END_VERSIONS
    """
}
