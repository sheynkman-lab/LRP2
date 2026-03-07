process NOVEL_PEPTIDES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(search_results)
    path novel_peptides_script
    path lr_cds_gtf
    path lr_orf_fasta

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
    echo "=========================================================================="
    echo "            NOVEL PEPTIDES CLASSIFICATION: ${meta.id}"
    echo "=========================================================================="
    echo "Sample: ${prefix}"
    echo "Search software: ${ms_search_software}"
    echo "Acquisition type: ${acquisition_type}"
    echo "LR CDS GTF: ${lr_cds_gtf}"
    echo "LR ORF FASTA: ${lr_orf_fasta}"
    echo "Input directory: \$(pwd)"
    echo "=========================================================================="
    echo ""

    # Copy search results to expected directory structure
    # The R script expects results in S5_PROTEOMICS/M2_FRAGPIPE/<sample_name>/ directory
    mkdir -p S5_PROTEOMICS/M2_FRAGPIPE/${prefix}

    # Copy all input files to the expected location
    if [ -d "${search_results}" ]; then
        cp -r ${search_results}/* S5_PROTEOMICS/M2_FRAGPIPE/${prefix}/ 2>/dev/null || true
    fi

    # Also check if files are directly provided (not in subdirectory)
    cp ${search_results}/*.tsv S5_PROTEOMICS/M2_FRAGPIPE/${prefix}/ 2>/dev/null || true

    echo "Directory structure created:"
    ls -lah S5_PROTEOMICS/M2_FRAGPIPE/${prefix}/ | head -20
    echo ""

    # Run novel peptides classification
    Rscript ${novel_peptides_script} \\
        --sample_name ${prefix} \\
        --ms_search_software ${ms_search_software} \\
        --acquisition_type ${acquisition_type} \\
        --lr_cds_gtf ${lr_cds_gtf} \\
        --lr_orf_fasta ${lr_orf_fasta} \\
        --outdir . \\
        $args

    # Verify output was created
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
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
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
