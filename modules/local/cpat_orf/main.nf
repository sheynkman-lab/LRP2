process CPAT_ORF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(corrected_fasta)
    path hexamer_file
    path logit_model

    output:
    tuple val(meta), path("*.predicted_proteome.CPAT.ORF_prob.tsv"), emit: orf_prob
    tuple val(meta), path("*.predicted_proteome.CPAT.ORF_prob.best.tsv"), emit: orf_prob_best
    tuple val(meta), path("*.predicted_proteome.CPAT.ORF_seqs.fa"), emit: orf_seqs
    tuple val(meta), path("*.predicted_proteome.CPAT.no_ORF.txt"), emit: no_orf
    tuple val(meta), path("*.predicted_proteome.CPAT.error"), emit: error_log
    tuple val(meta), path("*_S3_PREDICTED_PROTEOME_M1_CPAT_ORF_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_orf = task.ext.min_orf ?: params.min_orf
    def top_orf = task.ext.top_orf ?: params.top_orf

    """
    exec > >(tee ${prefix}_S3_PREDICTED_PROTEOME_M1_CPAT_ORF_log.txt) 2>&1

    cpat.py \\
        -x $hexamer_file \\
        -d $logit_model \\
        -g $corrected_fasta \\
        --min-orf=$min_orf \\
        --top-orf=$top_orf \\
        -o ${prefix}.predicted_proteome.CPAT \\
        $args \\
        2> ${prefix}.predicted_proteome.CPAT.error

    # Clean up CPAT run info log if it exists
    rm -f CPAT_run_info.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpat: \$(cpat.py --version 2>&1 | grep -oP 'CPAT-\\K[0-9.]+' || echo "3.0.4")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.predicted_proteome.CPAT.ORF_prob.tsv
    touch ${prefix}.predicted_proteome.CPAT.ORF_prob.best.tsv
    touch ${prefix}.predicted_proteome.CPAT.ORF_seqs.fa
    touch ${prefix}.predicted_proteome.CPAT.no_ORF.txt
    touch ${prefix}.predicted_proteome.CPAT.error
    touch ${prefix}_S3_PREDICTED_PROTEOME_M1_CPAT_ORF_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpat: 3.0.4
    END_VERSIONS
    """
}
