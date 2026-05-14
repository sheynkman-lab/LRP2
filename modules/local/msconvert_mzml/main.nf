process MSCONVERT_MZML {
    tag "$meta.id"
    label 'process_medium'
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:latest"

    input:
    tuple val(meta), path(raw_file)

    output:
    tuple val(meta), path("*.mzML"), emit: mzml
    tuple val(meta), path("*_S5_PROTEOMICS_M2_MSCONVERT_MZML_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def peak_picking = task.ext.peak_picking ?: params.msconvert_peak_picking ?: true
    def peak_picking_filter = peak_picking ? '--filter "peakPicking true 1-"' : ''
    def output_name = "${prefix}.mzML"

    """
    exec > >(tee ${prefix}_S5_PROTEOMICS_M2_MSCONVERT_MZML_log.txt) 2>&1

    wine msconvert \\
        $raw_file \\
        --outdir . \\
        $peak_picking_filter \\
        $args

    # Rename output to include prefix (sample name)
    # msconvert uses the input filename, so we need to rename
    original_output=\$(ls *.mzML | head -1)
    if [ "\$original_output" != "${output_name}" ]; then
        mv "\$original_output" "${output_name}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msconvert: \$(msconvert 2>&1 | grep -oP 'ProteoWizard release: \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mzML
    touch ${prefix}_S5_PROTEOMICS_M2_MSCONVERT_MZML_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msconvert: 3.0.0
    END_VERSIONS
    """
}
