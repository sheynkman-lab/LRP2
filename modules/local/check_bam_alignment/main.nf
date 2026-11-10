process CHECK_BAM_ALIGNMENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), stdout, emit: check
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Check if BAM header contains @PG line with pbmm2
    if samtools view -H ${bam} | grep -q '@PG.*pbmm2'; then
        echo "Sample '${meta.id}' will skip ISOCALL_ALIGN, as it is already a pbmm2-aligned BAM." >&2
        echo -n "true"
    else
        echo -n "false"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | sed -n '1s/samtools //p')
    END_VERSIONS
    """

    stub:
    """
    echo -n "false"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | sed -n '1s/samtools //p')
    END_VERSIONS
    """
}
