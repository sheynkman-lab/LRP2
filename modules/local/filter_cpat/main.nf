process FILTER_CPAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(orf_prob), path(orf_seqs), path(corrected_fasta), path(corrected_gtf), path(mapping_file)
    path reference_gtf
    path filter_cpat_script

    output:
    tuple val(meta), path("*_predicted_proteome_all_orfs_mapped.tsv"), emit: all_orfs_mapped
    tuple val(meta), path("*_predicted_proteome_corrected_filtered_CDS.gtf"), emit: cds_gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cpat_coding_threshold = task.ext.cpat_coding_threshold ?: params.cpat_coding_threshold

    """
    # Create output directory
    mkdir -p orf_calling

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    # Run the CPAT filter script with CLI arguments
    Rscript ${filter_cpat_script} \\
        --basename ${prefix} \\
        --cpat_fasta ${orf_seqs} \\
        --cpat_results ${orf_prob} \\
        --sample_fasta ${corrected_fasta} \\
        --sample_gtf ${corrected_gtf} \\
        --reference_gtf ${reference_gtf} \\
        --mapping_file ${mapping_file} \\
        --output_dir orf_calling \\
        --coding_threshold ${cpat_coding_threshold} \\
        $args

    # Move and rename outputs from orf_calling directory to match expected naming convention
    mv orf_calling/${prefix}.predicted_proteome.all_orfs_mapped.tsv ${prefix}_predicted_proteome_all_orfs_mapped.tsv
    mv orf_calling/${prefix}.predicted_proteome.corrected_filtered_CDS.gtf ${prefix}_predicted_proteome_corrected_filtered_CDS.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_proteome_all_orfs_mapped.tsv
    touch ${prefix}_predicted_proteome_corrected_filtered_CDS.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}
