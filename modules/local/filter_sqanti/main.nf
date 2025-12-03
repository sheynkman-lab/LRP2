process FILTER_SQANTI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:1.2"

    input:
    tuple val(meta), path(classification_file), path(corrected_gtf), path(corrected_fasta)
    path reference_gtf
    path reference_fasta
    path sample_metadata
    path filter_script
    path hashlib_script

    output:
    tuple val(meta), path("*_sqanti_transcript_classification_filtered.txt"), emit: classification_filtered
    tuple val(meta), path("*_sqanti_transcript_corrected_filtered.gtf"), emit: corrected_gtf_filtered
    tuple val(meta), path("*_sqanti_transcript_corrected_filtered.fasta"), emit: corrected_fasta_filtered
    tuple val(meta), path("*_sqanti_transcript_hashids_with_cpm_filtered.txt"), emit: hashids_filtered
    tuple val(meta), path("*_sqanti_transcript_all_hashids_with_cpm.txt"), emit: hashids_all
    tuple val(meta), path("*_sqanti_transcript_dropout_transcripts.tsv"), emit: dropout_transcripts
    tuple val(meta), path("*_sqanti_transcript_corrected_dropout.fasta"), emit: corrected_dropout_fasta
    tuple val(meta), path("*_sqanti_transcript_corrected_dropout.gtf"), emit: corrected_dropout_gtf
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
    def structure_filter = task.ext.structure_filter ?: params.structure_filter

    """
    # Link the SQANTI_QC output files to the working directory with expected naming
    ln -s \$(pwd)/$classification_file ${prefix}_sqanti_transcript_classification.txt
    ln -s \$(pwd)/$corrected_gtf ${prefix}_sqanti_transcript_corrected.gtf
    ln -s \$(pwd)/$corrected_fasta ${prefix}_sqanti_transcript_corrected.fasta

    # Export environment variables expected by the R script
    # OUTPUT_DIR must be the current working directory where symlinks and outputs are
    export OUTPUT_DIR=.
    export OUTPUT_BASE_NAME=$prefix
    export GENCODE_GTF_FILE=\$(pwd)/$reference_gtf
    export SAMPLE_METADATA=\$(pwd)/$sample_metadata
    export HASHLIB_SCRIPT=\$(pwd)/$hashlib_script
    export PROTEIN_CODING_FILTER=$protein_coding_filter
    export INTERNAL_PRIMING_FILTER=$internal_priming_filter
    export TEMPLATE_SWITCHING_FILTER=$template_switching_filter
    export PERCENT_POLYA_THRESHOLD=$percent_polya_threshold
    export STRUCTURE_FILTER=$structure_filter

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript \$(pwd)/$filter_script $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sqanti_transcript_classification_filtered.txt
    touch ${prefix}_sqanti_transcript_corrected_filtered.gtf
    touch ${prefix}_sqanti_transcript_corrected_filtered.fasta
    touch ${prefix}_sqanti_transcript_hashids_with_cpm_filtered.txt
    touch ${prefix}_sqanti_transcript_all_hashids_with_cpm.txt
    touch ${prefix}_sqanti_transcript_dropout_transcripts.tsv
    touch ${prefix}_sqanti_transcript_corrected_dropout.fasta
    touch ${prefix}_sqanti_transcript_corrected_dropout.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}
