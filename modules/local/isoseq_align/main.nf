process ISOSEQ_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.14.99--h9ee0642_0' :
        'quay.io/biocontainers/pbmm2:1.14.99--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    path(reference_fasta)

    output:
    tuple val(meta), path("*.aligned.bam"), emit: bam
    tuple val(meta), path("*.aligned.bam.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--preset ISOSEQ --sort'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbmm2 align \\
        $args \\
        $reference_fasta \\
        $bam \\
        ${prefix}.aligned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$( pbmm2 --version 2>&1 | sed 's/pbmm2 //;q' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aligned.bam
    touch ${prefix}.aligned.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$( pbmm2 --version 2>&1 | sed 's/pbmm2 //;q' )
    END_VERSIONS
    """
}

