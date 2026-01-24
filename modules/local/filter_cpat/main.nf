process FILTER_CPAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(orf_prob), path(orf_seqs), path(corrected_fasta), path(corrected_gtf)
    path reference_gtf
    path filter_cpat_script

    output:
    tuple val(meta), path("*_predicted_proteome_all_orfs_mapped.tsv"), emit: all_orfs_mapped
    tuple val(meta), path("*_predicted_proteome_best_orfs_mapped.tsv"), emit: best_orfs_mapped
    tuple val(meta), path("*_predicted_proteome_corrected_filtered_CDS.gtf"), emit: cds_gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cpat_coding_threshold = task.ext.cpat_coding_threshold ?: params.cpat_coding_threshold

    """
    # Create directory structure expected by the R script
    mkdir -p orf_calling sqanti_transcript

    # Create symlinks in the appropriate directories
    ln -sf \$(pwd)/$orf_prob orf_calling/${prefix}_cpat.ORF_prob.tsv
    ln -sf \$(pwd)/$orf_seqs orf_calling/${prefix}_cpat.ORF_seqs.fa
    ln -sf \$(pwd)/$corrected_fasta sqanti_transcript/${prefix}_corrected_filtered.fasta
    ln -sf \$(pwd)/$corrected_gtf sqanti_transcript/${prefix}_corrected_filtered.gtf

    # Export environment variables expected by the R script
    export OUTPUT_DIR=\$(pwd)
    export OUTPUT_BASE_NAME=$prefix
    export GENCODE_GTF_FILE=\$(pwd)/$reference_gtf
    export CPAT_CODING_THRESHOLD=$cpat_coding_threshold

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript \$(pwd)/$filter_cpat_script $args

    # Move and rename outputs from orf_calling directory to match expected naming convention
    mv orf_calling/${prefix}_all_orfs_mapped.tsv ${prefix}_predicted_proteome_all_orfs_mapped.tsv
    mv orf_calling/${prefix}_best_orfs_mapped.tsv ${prefix}_predicted_proteome_best_orfs_mapped.tsv
    mv orf_calling/${prefix}_corrected_filtered_CDS.gtf ${prefix}_predicted_proteome_corrected_filtered_CDS.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_proteome_all_orfs_mapped.tsv
    touch ${prefix}_predicted_proteome_best_orfs_mapped.tsv
    touch ${prefix}_predicted_proteome_corrected_filtered_CDS.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}
