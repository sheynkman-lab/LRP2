process LEAFCUTTER_LONGREAD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(filtered_gtf), path(transcript_counts)
    path sample_metadata
    path lr_leafcutter_script
    val control_group
    val min_samples_per_intron
    val min_samples_per_group
    val min_usage_ratio

    output:
    tuple val(meta), path("*_transcriptome_corrected_filtered.psl"), emit: psl
    tuple val(meta), path("*_transcriptome_corrected_filtered_intron_coords.txt"), emit: intron_coords
    tuple val(meta), path("*_transcriptome_corrected_filtered_exon_coords.txt"), emit: exon_coords
    tuple val(meta), path("*_lr_leafcutter_subisoform_clusters.txt"), emit: subisoform_clusters
    tuple val(meta), path("*_lr_leafcutter_perind_numers.counts.gz"), emit: counts_matrix
    tuple val(meta), path("*_groups_file.txt"), emit: groups_file
    tuple val(meta), path("Rplots.pdf"), emit: rplots, optional: true
    tuple val(meta), path("*_cluster_significance.txt"), emit: cluster_significance, optional: true
    tuple val(meta), path("*_effect_sizes.txt"), emit: effect_sizes, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_samp_intron = min_samples_per_intron ?: 2
    def min_samp_group = min_samples_per_group ?: 2
    def usage_ratio = min_usage_ratio ?: 0.01

    """
    mkdir -p multisample_analysis/lr_leafcutter

    export OUTPUT_BASE_NAME=${prefix}
    export OUTPUT_DIR=.
    export SAMPLE_METADATA=\$(pwd)/${sample_metadata}

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    Rscript ${lr_leafcutter_script} \\
        --gtf ${filtered_gtf} \\
        --counts ${transcript_counts} \\
        --mode exon

    # Move leafcutter outputs from subdirectory to work directory root so detected
    mv multisample_analysis/lr_leafcutter/* . 2>/dev/null || true

    # Create Python virtual environment for leafcutter differential splicing
    python3 -m venv leafcutter_env
    source leafcutter_env/bin/activate

    # Install leafcutter requirements
    pip install --upgrade pip
    pip install numpy scipy pandas rpy2

    # Run differential splicing if leafcutter_ds.py script is available
    if [ -f "\${LEAFCUTTER_DS_SCRIPT:-/dev/null}" ]; then
        python \${LEAFCUTTER_DS_SCRIPT} \\
            -0 ${control_group} \\
            --min_samples_per_intron ${min_samp_intron} \\
            --min_samples_per_group ${min_samp_group} \\
            --num_threads ${task.cpus} \\
            -o ${prefix} \\
            ${prefix}_lr_leafcutter_perind_numers.counts.gz \\
            ${prefix}_groups_file.txt || echo "Leafcutter DS step skipped or failed"
    fi

    deactivate

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_transcriptome_corrected_filtered.psl
    touch ${prefix}_transcriptome_corrected_filtered_intron_coords.txt
    touch ${prefix}_transcriptome_corrected_filtered_exon_coords.txt
    touch ${prefix}_lr_leafcutter_subisoform_clusters.txt
    touch ${prefix}_lr_leafcutter_perind_numers.counts.gz
    touch ${prefix}_groups_file.txt
    touch Rplots.pdf
    touch ${prefix}_cluster_significance.txt
    touch ${prefix}_effect_sizes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.5.0
        python: 3.11.4
    END_VERSIONS
    """
}
