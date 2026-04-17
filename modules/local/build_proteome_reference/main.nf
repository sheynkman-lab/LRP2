process BUILD_PROTEOME_REFERENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(lrp_fasta), path(counts), path(custom_fasta), path(gencode_protein_fasta)
    path build_proteome_reference_script
    val genome_name
    val no_gencode

    output:
    tuple val(meta), path("*.proteomics.reference.fasta"), emit: reference_fasta
    tuple val(meta), path("*.proteomics.reference.tsv"), emit: reference_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Handle both generic NO_FILE and sample-specific placeholder names (e.g., sample1_NO_LRP_FASTA)
    def lrp_fasta_arg = (lrp_fasta.name != 'NO_FILE' && !lrp_fasta.name.contains('_NO_LRP_FASTA')) ? "--lrp_fasta ${lrp_fasta}" : ""
    def counts_arg = (counts.name != 'NO_FILE' && !counts.name.contains('_NO_COUNTS')) ? "--counts ${counts}" : ""
    //def gencode_fasta_arg = "--gencode_fasta ${gencode_fasta}"
    def gencode_fasta_arg = (gencode_protein_fasta.name != 'NO_FILE' && !gencode_protein_fasta.name.contains('_NO_GENCODE_PROTEIN_FASTA')) ? "--gencode_fasta ${gencode_protein_fasta}" : ""
    def gencode_flag = no_gencode ? "--no_gencode" : ""
    def custom_fasta_arg = (custom_fasta.name != 'NO_FILE' && !custom_fasta.name.contains('_NO_CUSTOM_FASTA')) ? "--custom_fasta ${custom_fasta}" : ""

    """
    Rscript ${build_proteome_reference_script} \\
        ${lrp_fasta_arg} \\
        ${counts_arg} \\
        ${custom_fasta_arg} \\
        ${gencode_fasta_arg} \\
        --genome_name ${genome_name} \\
        ${gencode_flag} \\
        --sample_name ${prefix} \\
        --outdir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        bioconductor-biostrings: \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch test.proteomics.reference.fasta
    touch test.proteomics.reference.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.5.0
        bioconductor-biostrings: 2.72.0
    END_VERSIONS
    """
}
