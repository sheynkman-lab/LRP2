process LEAFCUTTER_LONGREAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(filtered_gtf), path(transcript_counts)
    path sample_metadata
    path lr_leafcutter_script
    path leafcutter_ds_script
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
    def prefix = task.ext.prefix ?: "${control_group}_${experimental_group}"
    def min_samp_intron = min_samples_per_intron ?: 2
    def min_samp_group = min_samples_per_group ?: 2
    def usage_ratio = min_usage_ratio ?: 0.01

    """
    mkdir -p multisample_analysis/lr_leafcutter

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    # Run leafcutter clustering
    Rscript ${lr_leafcutter_script} \\
        --gtf ${filtered_gtf} \\
        --counts ${transcript_counts} \\
        --sample_metadata ${sample_metadata} \\
        --basename ${prefix} \\
        --output_dir multisample_analysis/lr_leafcutter \\
        --mode exon \\
        --min_usage_ratio ${usage_ratio} \\
        $args

    # Move leafcutter outputs from subdirectory to work directory root so detected
    mv multisample_analysis/lr_leafcutter/* . 2>/dev/null || true

    # Run differential splicing if script is provided
    # Note: leafcutter_ds.py requires the 'leafcutter' Python package which should be
    # available in the PYTHONPATH from the bin directory
    if [ -f "${leafcutter_ds_script}" ]; then
        # Add bin directory to Python path so leafcutter package can be imported
        export PYTHONPATH=\$(dirname ${leafcutter_ds_script}):\${PYTHONPATH:-}

        python ${leafcutter_ds_script} \\
            -0 ${control_group} \\
            --min_samples_per_intron ${min_samp_intron} \\
            --min_samples_per_group ${min_samp_group} \\
            --num_threads ${task.cpus} \\
            -o ${prefix} \\
            ${prefix}_lr_leafcutter_perind_numers.counts.gz \\
            ${prefix}_groups_file.txt || echo "WARNING: Leafcutter differential splicing failed"
    else
        echo "NOTE: Leafcutter differential splicing script not provided - skipping DS analysis"
        echo "      Clustering outputs (PSL, coords, subisoform counts) are still available"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        python: \$(python3 --version | sed 's/Python //g')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${control_group}_${experimental_group}"
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
