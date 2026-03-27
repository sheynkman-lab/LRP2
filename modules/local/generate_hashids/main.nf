process GENERATE_HASHIDS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(corrected_gtf), path(classification)
    path reference_gtf
    path hashlib_script
    path generate_hashids_script

    output:
    tuple val(meta), path("*.transcriptome.hashids_mapping.txt"), emit: hashids_mapping
    tuple val(meta), path("*.transcriptome.psl"), emit: corrected_psl
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    
    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript ${generate_hashids_script} \\
        --basename ${prefix} \\
        --sample_gtf ${corrected_gtf} \\
        --classification ${classification} \\
        --reference_gtf ${reference_gtf} \\
        --hashlib_script ${hashlib_script} \\
        --output_dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcriptome.hashids_mapping.txt
    touch ${prefix}.transcriptome.psl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
        python: 3.11.0
    END_VERSIONS
    """
}
