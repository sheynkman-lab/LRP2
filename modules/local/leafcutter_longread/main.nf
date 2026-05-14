process LEAFCUTTER_LONGREAD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/leafcutter-longread:latest"

    input:
    tuple val(meta), path(filtered_gtf), path(transcript_counts)
    path sample_metadata
    path lr_leafcutter_script
    path leafcutter_ds_script
    val control_group
    val experimental_group
    val min_samples_per_intron
    val min_samples_per_group
    val min_usage_ratio
    val dataset_name
    val leafcutter_threads

    output:
    tuple val(meta), path("*transcriptome.filtered.psl"), emit: psl
    tuple val(meta), path("*transcriptome.filtered_intron_coords.txt"), emit: intron_coords
    tuple val(meta), path("*transcriptome.filtered_exon_coords.txt"), emit: exon_coords
    tuple val(meta), path("*.lr_leafcutter.subisoform_clusters.txt"), emit: subisoform_clusters
    tuple val(meta), path("*.lr_leafcutter.perind_numers.counts.gz"), emit: counts_matrix
    tuple val(meta), path("*.lr_leafcutter.groups_file.txt"), emit: groups_file
    tuple val(meta), path("lr_leafcutter.Rplots.pdf"), emit: rplots, optional: true
    tuple val(meta), path("*.lr_leafcutter.ds_cluster_significance.txt"), emit: cluster_significance, optional: true
    tuple val(meta), path("*.lr_leafcutter.ds_effect_sizes.txt"), emit: effect_sizes, optional: true
    tuple val(meta), path("*_S4_MULTISAMPLE_ANALYSIS_M1_LEAFCUTTER_LONGREAD_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_samp_intron = min_samples_per_intron ?: 2
    def min_samp_group = min_samples_per_group ?: 2
    def usage_ratio = min_usage_ratio ?: 0.01
    def threads = leafcutter_threads ?: task.cpus

    """
    exec > >(tee ${prefix}_S4_MULTISAMPLE_ANALYSIS_M1_LEAFCUTTER_LONGREAD_log.txt) 2>&1

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    # Run leafcutter clustering
    Rscript ${lr_leafcutter_script} \\
        --gtf ${filtered_gtf} \\
        --counts ${transcript_counts} \\
        --sample_metadata ${sample_metadata} \\
        --basename ${prefix} \\
        --output_dir . \\
        --mode exon \\
        --min_usage_ratio ${usage_ratio} \\
        $args

    # Run differential splicing if script is provided
    # Note: leafcutter_ds.py requires the 'leafcutter' Python package which should be available in the PYTHONPATH from the bin directory
    if [ -f "${leafcutter_ds_script}" ]; then
    
        # Filter groups file to only include samples from the two groups being compared since leafcutter_ds.py only supports pairwise comparisons
        awk -v ctrl="${control_group}" -v exper="${experimental_group}" '\$2 == ctrl || \$2 == exper' ${prefix}.lr_leafcutter.groups_file.txt > ${prefix}.lr_leafcutter.filtered_groups_file.txt
        n_control=\$(awk -v ctrl="${control_group}" '\$2 == ctrl' ${prefix}.lr_leafcutter.filtered_groups_file.txt | wc -l)
        n_experimental=\$(awk -v exper="${experimental_group}" '\$2 == exper' ${prefix}.lr_leafcutter.filtered_groups_file.txt | wc -l)

        if [ \$n_control -ge ${min_samp_group} ] && [ \$n_experimental -ge ${min_samp_group} ]; then

            export PYTHONPATH=\$(dirname ${leafcutter_ds_script}):\${PYTHONPATH:-}
            python ${leafcutter_ds_script} \\
                -0 ${control_group} \\
                --min_samples_per_intron ${min_samp_intron} \\
                --min_samples_per_group ${min_samp_group} \\
                --num_threads ${threads} \\
                -o ${prefix} \\
                ${prefix}.lr_leafcutter.perind_numers.counts.gz \\
                ${prefix}.lr_leafcutter.filtered_groups_file.txt || echo "WARNING: Leafcutter differential splicing failed"

            if [ -f "${prefix}_cluster_significance.txt" ]; then
                mv ${prefix}_cluster_significance.txt ${prefix}.lr_leafcutter.ds_cluster_significance.txt
            fi
            if [ -f "${prefix}_effect_sizes.txt" ]; then
                mv ${prefix}_effect_sizes.txt ${prefix}.lr_leafcutter.ds_effect_sizes.txt
            fi

            # Replace the unfiltered groups file with the filtered one used by analysis for output dir
            mv ${prefix}.lr_leafcutter.filtered_groups_file.txt ${prefix}.lr_leafcutter.groups_file.txt

        else
            echo "WARNING: Insufficient number of samples for differential splicing analysis! "
            echo "         leafcutter-longread requires at least ${min_samp_group} samples per group."
            echo "         Found only \$n_control ${control_group} and \$n_experimental ${experimental_group} samples." 
            echo "         Please adjust your sample groups or reduce the minimum samples per group threshold to run differential splicing analysis."
        fi
    else
        echo "NOTE: Leafcutter differential splicing script not provided - skipping DS analysis."
        echo "      Clustering outputs (PSL, coords, subisoform counts) will still be generated."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${dataset_name}.transcriptome.filtered.psl
    touch ${dataset_name}.transcriptome.filtered_intron_coords.txt
    touch ${dataset_name}.transcriptome.filtered_exon_coords.txt
    touch ${prefix}.lr_leafcutter.subisoform_clusters.txt
    touch ${prefix}.lr_leafcutter.perind_numers.counts.gz
    touch ${prefix}.lr_leafcutter.groups_file.txt
    touch lr_leafcutter.Rplots.pdf
    touch ${prefix}.lr_leafcutter.ds_cluster_significance.txt
    touch ${prefix}.lr_leafcutter.ds_effect_sizes.txt
    touch ${prefix}_S4_MULTISAMPLE_ANALYSIS_M1_LEAFCUTTER_LONGREAD_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.5.0
        python: 3.11.4
    END_VERSIONS
    """
}
