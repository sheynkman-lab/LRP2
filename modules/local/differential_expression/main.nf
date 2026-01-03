process DIFFERENTIAL_EXPRESSION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(transcript_counts), path(orf_counts)
    path sample_metadata
    path multisample_script
    val control_group
    val experimental_group

    output:
    tuple val(meta), path("differential_gene_expression/*_DGE_edgeR_results.txt"), emit: dge_results
    tuple val(meta), path("differential_gene_expression/*_DGE_edgeR_raw_CPM_matrix.txt"), emit: dge_raw_cpm
    tuple val(meta), path("differential_gene_expression/*_DGE_edgeR_normalized_CPM_matrix.txt"), emit: dge_norm_cpm
    tuple val(meta), path("differential_gene_expression/*_DGE_MD_plot.pdf"), emit: dge_plot
    tuple val(meta), path("differential_transcript_expression/*_DTE_edgeR_results.txt"), emit: dte_results
    tuple val(meta), path("differential_transcript_expression/*_DTE_edgeR_raw_CPM_matrix.txt"), emit: dte_raw_cpm
    tuple val(meta), path("differential_transcript_expression/*_DTE_edgeR_normalized_CPM_matrix.txt"), emit: dte_norm_cpm
    tuple val(meta), path("differential_transcript_expression/*_DTE_MD_plot.pdf"), emit: dte_plot
    tuple val(meta), path("differential_ORF_expression/*_DE_ORF_edgeR_results.txt"), emit: de_orf_results
    tuple val(meta), path("differential_ORF_expression/*_DE_ORF_edgeR_raw_CPM_matrix.txt"), emit: de_orf_raw_cpm
    tuple val(meta), path("differential_ORF_expression/*_DE_ORF_edgeR_normalized_CPM_matrix.txt"), emit: de_orf_norm_cpm
    tuple val(meta), path("differential_ORF_expression/*_DE_ORF_MD_plot.pdf"), emit: de_orf_plot
    tuple val(meta), path("differential_transcript_usage/*_DTU_transcript_DRIMSeq_summary.txt"), emit: dtu_summary
    tuple val(meta), path("differential_ORF_usage/*_DU_ORF_DRIMSeq_summary.txt"), emit: du_orf_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Output subdirectories expected by R script
    mkdir -p differential_gene_expression
    mkdir -p differential_transcript_expression
    mkdir -p differential_ORF_expression
    mkdir -p differential_transcript_usage
    mkdir -p differential_ORF_usage
    mkdir -p sqanti_transcript
    mkdir -p protein_sqanti
    mkdir -p multisample_analysis

    export OUTPUT_BASE_NAME=${prefix}
    export OUTPUT_DIR=.
    export SAMPLE_METADATA=\$(pwd)/${sample_metadata}
    export CONTROL_GROUP=${control_group}
    export EXPERIMENTAL_GROUP=${experimental_group}

    # Link input files to expected names and locations for the R script
    ln -sf \$(pwd)/${transcript_counts} sqanti_transcript/${prefix}_hashids_with_cpm_filtered.txt
    ln -sf \$(pwd)/${orf_counts} protein_sqanti/${prefix}_hashids_with_cpm_ORF.txt

    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript ${multisample_script} $args

    # Move outputs from multisample_analysis to organized subdirectories
    mv multisample_analysis/${prefix}_DGE_* differential_gene_expression/ 2>/dev/null || true
    mv multisample_analysis/${prefix}_DTE_* differential_transcript_expression/ 2>/dev/null || true
    mv multisample_analysis/${prefix}_DE_ORF_* differential_ORF_expression/ 2>/dev/null || true
    mv multisample_analysis/${prefix}_DTU_* differential_transcript_usage/ 2>/dev/null || true
    mv multisample_analysis/${prefix}_DU_ORF_* differential_ORF_usage/ 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        bioconductor-edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
        bioconductor-drimseq: \$(Rscript -e "cat(as.character(packageVersion('DRIMSeq')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p differential_gene_expression
    mkdir -p differential_transcript_expression
    mkdir -p differential_ORF_expression
    mkdir -p differential_transcript_usage
    mkdir -p differential_ORF_usage

    touch differential_gene_expression/${prefix}_DGE_edgeR_results.txt
    touch differential_gene_expression/${prefix}_DGE_edgeR_raw_CPM_matrix.txt
    touch differential_gene_expression/${prefix}_DGE_edgeR_normalized_CPM_matrix.txt
    touch differential_gene_expression/${prefix}_DGE_MD_plot.pdf

    touch differential_transcript_expression/${prefix}_DTE_edgeR_results.txt
    touch differential_transcript_expression/${prefix}_DTE_edgeR_raw_CPM_matrix.txt
    touch differential_transcript_expression/${prefix}_DTE_edgeR_normalized_CPM_matrix.txt
    touch differential_transcript_expression/${prefix}_DTE_MD_plot.pdf

    touch differential_ORF_expression/${prefix}_DE_ORF_edgeR_results.txt
    touch differential_ORF_expression/${prefix}_DE_ORF_edgeR_raw_CPM_matrix.txt
    touch differential_ORF_expression/${prefix}_DE_ORF_edgeR_normalized_CPM_matrix.txt
    touch differential_ORF_expression/${prefix}_DE_ORF_MD_plot.pdf

    touch differential_transcript_usage/${prefix}_DTU_transcript_DRIMSeq_summary.txt
    touch differential_ORF_usage/${prefix}_DU_ORF_DRIMSeq_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.5.0
        bioconductor-edger: 4.2.0
        bioconductor-drimseq: 1.32.0
    END_VERSIONS
    """
}
