process PROTEIN_CLASSIFICATION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(protein_classification), path(cds_gtf), path(corrected_fasta), path(hashids_cpm), path(all_orfs_mapped)
    path reference_gtf
    path protein_class_script

    output:
    tuple val(meta), path("*.predicted_proteome.best_ORF_summary.txt"), emit: protein_all_isoforms
    tuple val(meta), path("*.predicted_proteome.best_ORF.fa"), emit: protein_all_orfs_fasta
    tuple val(meta), path("*.predicted_proteome.collapsed_high_confidence_ORF_hashids_with_cpm.txt"), emit: hashids_orf
    tuple val(meta), path("*.predicted_proteome.collapsed_high_confidence_ORF.gtf"), emit: protein_gtf
    tuple val(meta), path("*.predicted_proteome.collapsed_high_confidence_ORF.bed"), emit: protein_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_junctions_after_stop = task.ext.min_junctions_after_stop ?: params.min_junctions_after_stop_codon
    def protein_class_keep = task.ext.protein_class_keep ?: params.protein_class_keep

    """
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript \$(pwd)/$protein_class_script \\
        --basename $prefix \\
        --gencode_gtf $reference_gtf \\
        --sample_cds_gtf $cds_gtf \\
        --sample_dna_fasta $corrected_fasta \\
        --mapped_orfs $all_orfs_mapped \\
        --protein_sqanti $protein_classification \\
        --cpm_file $hashids_cpm \\
        --output_dir . \\
        --min_junctions_after_stop $min_junctions_after_stop \\
        --protein_class_keep "$protein_class_keep" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.predicted_proteome.best_ORF_summary.txt
    touch ${prefix}.predicted_proteome.best_ORF.fa
    touch ${prefix}.predicted_proteome.collapsed_high_confidence_ORF_hashids_with_cpm.txt
    touch ${prefix}.predicted_proteome.collapsed_high_confidence_ORF.gtf
    touch ${prefix}.predicted_proteome.collapsed_high_confidence_ORF.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}
