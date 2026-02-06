process FILTER_TRANSCRIPTOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(classification_file), path(corrected_gtf), path(corrected_fasta), path(hashids_mapping), path(corrected_psl)
    path reference_gtf
    path reference_fasta
    path sample_metadata
    path filter_script

    output:
    tuple val(meta), path("*_transcriptome_classification.txt"), emit: classification
    tuple val(meta), path("*_transcriptome_classification_filtered.txt"), emit: classification_filtered
    tuple val(meta), path("*_transcriptome_corrected.gtf"), emit: corrected_gtf
    tuple val(meta), path("*_transcriptome_corrected_filtered.gtf"), emit: corrected_gtf_filtered
    tuple val(meta), path("*_transcriptome_corrected.fasta"), emit: corrected_fasta
    tuple val(meta), path("*_transcriptome_corrected_filtered.fasta"), emit: corrected_fasta_filtered
    tuple val(meta), path("*_transcriptome_corrected.psl"), emit: corrected_psl
    tuple val(meta), path("*_transcriptome_hashids_mapping.txt"), emit: hashids_mapping
    tuple val(meta), path("*_transcriptome_hashids_with_cpm_filtered.txt"), emit: hashids_filtered
    tuple val(meta), path("*_transcriptome_all_hashids_with_cpm.txt"), emit: hashids_all
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
    # Create directory structure expected by the R script
    mkdir -p sqanti_transcript

    # Link the SQANTI_QC output files and hashids_mapping with expected naming
    ln -sf \$(pwd)/$classification_file sqanti_transcript/${prefix}_classification.txt
    ln -sf \$(pwd)/$corrected_gtf sqanti_transcript/${prefix}_corrected.gtf
    ln -sf \$(pwd)/$corrected_fasta sqanti_transcript/${prefix}_corrected.fasta
    ln -sf \$(pwd)/$hashids_mapping sqanti_transcript/${prefix}_hashids_mapping.txt
    ln -sf \$(pwd)/$corrected_psl sqanti_transcript/${prefix}_corrected.psl

    # Export environment variables expected by the R script
    export OUTPUT_DIR=\$(pwd)
    export OUTPUT_BASE_NAME=$prefix
    export GENCODE_GTF_FILE=\$(pwd)/$reference_gtf
    export SAMPLE_METADATA=\$(pwd)/$sample_metadata
    export PROTEIN_CODING_FILTER=$protein_coding_filter
    export INTERNAL_PRIMING_FILTER=$internal_priming_filter
    export TEMPLATE_SWITCHING_FILTER=$template_switching_filter
    export PERCENT_POLYA_THRESHOLD=$percent_polya_threshold
    export TRANSCRIPT_CLASS_KEEP=$transcript_class_keep

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    # Run the filter transcriptome script (expects hashids_mapping.txt to already exist)
    Rscript \$(pwd)/$filter_script $args

    # Move outputs to expected locations with _transcriptome prefix
    mv sqanti_transcript/${prefix}_classification.txt ${prefix}_transcriptome_classification.txt
    mv sqanti_transcript/${prefix}_corrected.gtf ${prefix}_transcriptome_corrected.gtf
    mv sqanti_transcript/${prefix}_corrected.fasta ${prefix}_transcriptome_corrected.fasta
    mv sqanti_transcript/${prefix}_corrected.psl ${prefix}_transcriptome_corrected.psl
    mv sqanti_transcript/${prefix}_hashids_mapping.txt ${prefix}_transcriptome_hashids_mapping.txt
    mv sqanti_transcript/${prefix}_classification_filtered.txt ${prefix}_transcriptome_classification_filtered.txt
    mv sqanti_transcript/${prefix}_corrected_filtered.gtf ${prefix}_transcriptome_corrected_filtered.gtf
    mv sqanti_transcript/${prefix}_corrected_filtered.fasta ${prefix}_transcriptome_corrected_filtered.fasta
    mv sqanti_transcript/${prefix}_hashids_with_cpm_filtered.txt ${prefix}_transcriptome_hashids_with_cpm_filtered.txt
    mv sqanti_transcript/${prefix}_all_hashids_with_cpm.txt ${prefix}_transcriptome_all_hashids_with_cpm.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_transcriptome_classification.txt
    touch ${prefix}_transcriptome_classification_filtered.txt
    touch ${prefix}_transcriptome_corrected.gtf
    touch ${prefix}_transcriptome_corrected_filtered.gtf
    touch ${prefix}_transcriptome_corrected.fasta
    touch ${prefix}_transcriptome_corrected_filtered.fasta
    touch ${prefix}_transcriptome_corrected.psl
    touch ${prefix}_transcriptome_hashids_mapping.txt
    touch ${prefix}_transcriptome_hashids_with_cpm_filtered.txt
    touch ${prefix}_transcriptome_all_hashids_with_cpm.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}
