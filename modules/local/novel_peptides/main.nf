process NOVEL_PEPTIDES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(fragpipe_peptide_file), path(reference_fasta), path(custom_gtf), path(custom_fasta), path(gencode_gtf), path(gencode_fasta)
    path novel_peptides_script
    val gencode_version

    output:
    tuple val(meta), path("*.proteomics.novel_peptides.tsv"), emit: novel_peptides
    tuple val(meta), path("*.proteomics.all_peptides.bed"), optional: true, emit: peptides_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ms_search_software = meta.protein_search ?: 'fragpipe'
    def acquisition_type = meta.mass_spec_type ?: 'DDA'
    def has_custom_data = (!custom_gtf.name.contains('_NO_GTF') && !custom_fasta.name.contains('_NO_FASTA'))
    def custom_gtf_arg = has_custom_data ? "--custom_gtf ${custom_gtf}" : ""
    def custom_fasta_arg = has_custom_data ? "--custom_fasta ${custom_fasta}" : ""
    def has_gencode_data = (!gencode_gtf.name.contains('_NO_GTF') && !gencode_fasta.name.contains('_NO_FASTA'))
    def gencode_gtf_arg = has_gencode_data ? "--gencode_gtf ${gencode_gtf}" : ""
    def gencode_fasta_arg = has_gencode_data ? "--gencode_fasta ${gencode_fasta}" : ""

    """
    echo "NOVEL PEPTIDES CLASSIFICATION: ${meta.id}"
    echo "Sample: ${prefix}"
    echo "Search software: ${ms_search_software}"
    echo "Acquisition type: ${acquisition_type}"
    echo "Custom GTF: ${has_custom_data ? custom_gtf : 'N/A (no custom data)'}"
    echo "Custom FASTA: ${has_custom_data ? custom_fasta : 'N/A (no custom data)'}"
    echo "GENCODE GTF: ${has_gencode_data ? gencode_gtf : 'N/A (will use fallback)'}"
    echo "GENCODE FASTA: ${has_gencode_data ? gencode_fasta : 'N/A (will use fallback)'}"
    echo "FragPipe peptide file: ${fragpipe_peptide_file}"
    echo ""

    # Derive FragPipe results directory from the peptide file path
    if [ ! -f "${fragpipe_peptide_file}" ]; then
        echo "ERROR: FragPipe peptide file not found: ${fragpipe_peptide_file}"
        echo "Please check that FragPipe completed successfully."
        exit 1
    fi

    FRAGPIPE_PEPTIDE_REALPATH=\$(realpath "${fragpipe_peptide_file}")
    FRAGPIPE_RESULTS_DIR=\$(dirname "\$FRAGPIPE_PEPTIDE_REALPATH")
    echo "FragPipe results directory: \$FRAGPIPE_RESULTS_DIR"
    echo ""
    if [ ! -d "\$FRAGPIPE_RESULTS_DIR" ]; then
        echo "ERROR: FragPipe results directory not found: \$FRAGPIPE_RESULTS_DIR"
        exit 1
    fi

    Rscript ${novel_peptides_script} \\
        --sample_name ${prefix} \\
        --ms_search_software ${ms_search_software} \\
        --acquisition_type ${acquisition_type} \\
        --fragpipe_results_dir \$FRAGPIPE_RESULTS_DIR \\
        ${custom_gtf_arg} \\
        ${custom_fasta_arg} \\
        ${gencode_gtf_arg} \\
        ${gencode_fasta_arg} \\
        --gencode_version ${gencode_version} \\
        --outdir . \\
        $args

    if [ ! -f "${prefix}.proteomics.novel_peptides.tsv" ]; then
        echo "ERROR: Novel peptides output file not created"
        exit 1
    fi

    echo ""
    echo "Novel peptides classification completed successfully!"
    echo "Output: ${prefix}.proteomics.novel_peptides.tsv"
    echo ""

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>/dev/null | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g' || echo "unknown")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.proteomics.novel_peptides.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.5.0
        r-tidyverse: 2.0.0
    END_VERSIONS
    """
}
