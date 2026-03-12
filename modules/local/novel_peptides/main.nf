process NOVEL_PEPTIDES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), val(fragpipe_results_dir)
    path novel_peptides_script
    path lr_cds_gtf
    path lr_orf_fasta
    val gencode_version

    output:
    tuple val(meta), path("*_novel_peptides.tsv"), emit: novel_peptides
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ms_search_software = meta.protein_search ?: 'fragpipe'
    def acquisition_type = meta.mass_spec_type ?: 'DDA'

    """
    echo "NOVEL PEPTIDES CLASSIFICATION: ${meta.id}"
    echo "Sample: ${prefix}"
    echo "Search software: ${ms_search_software}"
    echo "Acquisition type: ${acquisition_type}"
    echo "LR CDS GTF: ${lr_cds_gtf}"
    echo "LR ORF FASTA: ${lr_orf_fasta}"
    echo "FragPipe results directory: ${fragpipe_results_dir}"
    echo ""

    if [ ! -d "${fragpipe_results_dir}" ]; then
        echo "ERROR: FragPipe results directory not found: ${fragpipe_results_dir}"
        echo "Please check that FragPipe completed successfully."
        exit 1
    fi

    Rscript ${novel_peptides_script} \\
        --sample_name ${prefix} \\
        --ms_search_software ${ms_search_software} \\
        --acquisition_type ${acquisition_type} \\
        --fragpipe_results_dir ${fragpipe_results_dir} \\
        --lr_cds_gtf ${lr_cds_gtf} \\
        --lr_orf_fasta ${lr_orf_fasta} \\
        --gencode_version ${gencode_version} \\
        --outdir . \\
        $args

    if [ ! -f "${prefix}_novel_peptides.tsv" ]; then
        echo "ERROR: Novel peptides output file not created"
        exit 1
    fi

    echo ""
    echo "Novel peptides classification completed successfully!"
    echo "Output: ${prefix}_novel_peptides.tsv"
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
    touch ${prefix}_novel_peptides.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.5.0
        r-tidyverse: 2.0.0
    END_VERSIONS
    """
}
