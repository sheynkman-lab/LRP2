process FILTER_TRANSCRIPTOME {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(classification_file), path(corrected_gtf), path(corrected_fasta), path(hashids_mapping), path(corrected_psl)
    path reference_gtf
    path reference_fasta
    path sample_metadata
    path filter_script

    output:
    tuple val(meta), path("*.transcriptome.filtered_SQANTI_classification.txt"), emit: classification_filtered
    tuple val(meta), path("*.transcriptome.filtered.gtf"), emit: corrected_gtf_filtered
    tuple val(meta), path("*.transcriptome.filtered.bed"), emit: corrected_bed_filtered
    tuple val(meta), path("*.transcriptome.filtered.fasta"), emit: corrected_fasta_filtered
    tuple val(meta), path("*.transcriptome.filtered_hashids_with_cpm.txt"), emit: hashids_filtered
    tuple val(meta), path("*.transcriptome.all_hashids_with_cpm.txt"), emit: hashids_all
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def protein_coding_filter = task.ext.protein_coding_filter ?: params.protein_coding_filter
    def internal_priming_filter = task.ext.internal_priming_filter ?: params.internal_priming_filter
    def template_switching_filter = task.ext.template_switching_filter ?: params.template_switching_filter
    def percent_polya_threshold = task.ext.percent_polya_threshold ?: params.percent_polya_threshold
    def transcript_class_keep = task.ext.transcript_class_keep ?: params.transcript_class_keep

    """
    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript ${filter_script} \\
        --basename ${prefix} \\
        --classification ${classification_file} \\
        --sample_gtf ${corrected_gtf} \\
        --sample_fasta ${corrected_fasta} \\
        --mapping_file ${hashids_mapping} \\
        --output_dir . \\
        --filter_protein_coding ${protein_coding_filter} \\
        --filter_internal_priming ${internal_priming_filter} \\
        --filter_rts ${template_switching_filter} \\
        --percent_polya_threshold ${percent_polya_threshold} \\
        --transcript_class_keep ${transcript_class_keep} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcriptome.filtered_SQANTI_classification.txt
    touch ${prefix}.transcriptome.filtered.gtf
    touch ${prefix}.transcriptome.filtered.bed
    touch ${prefix}.transcriptome.filtered.fasta
    touch ${prefix}.transcriptome.filtered_hashids_with_cpm.txt
    touch ${prefix}.transcriptome.all_hashids_with_cpm.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}
