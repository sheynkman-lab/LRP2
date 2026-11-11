process MSCONVERT_MZML {
    tag "$meta.id"
    label 'process_medium'
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "docker.io/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:latest"

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
    // Conda only: options for ThermoRawFileParser (ext.args holds msconvert options, which it doesn't accept)
    def args2 = task.ext.args2 ?: ''
    def trfp_peak_picking = peak_picking ? '' : '--noPeakPicking'

    """
    exec > >(tee ${prefix}_S5_PROTEOMICS_M2_MSCONVERT_MZML_log.txt) 2>&1

    # Docker/Singularity: msconvert from the ProteoWizard image (any vendor format), run under wine.
    # Conda: msconvert's vendor libraries are Windows-only and have no Conda package, so this module's
    # environment.yml provides ThermoRawFileParser instead. It reads Thermo .raw files only. The
    # ProteoWizard image doesn't contain it, so container runs always take the msconvert branch.
    if command -v ThermoRawFileParser > /dev/null 2>&1; then
        case "$raw_file" in
            *.raw|*.RAW) ;;
            *)
                echo "ERROR: with Conda, MSCONVERT_MZML converts Thermo .raw files only (ThermoRawFileParser)." >&2
                echo "       '$raw_file' is another format: run with Docker or Singularity, or supply mzML files." >&2
                exit 1
                ;;
        esac
        echo "Converting $raw_file with ThermoRawFileParser (Conda)"
        # --format=2: indexed mzML, msconvert's default output. Peaks are centroided with the
        # vendor algorithm unless --noPeakPicking is given, like msconvert's "peakPicking true 1-".
        ThermoRawFileParser \\
            --input=$raw_file \\
            --output=${output_name} \\
            --format=2 \\
            $trfp_peak_picking \\
            $args2
        CONVERTER_NAME="thermorawfileparser"
        CONVERTER_VERSION=\$(ThermoRawFileParser --version 2>&1 | grep -oE '[0-9]+(\\.[0-9]+)+' | head -n 1 || true)
    else
        wine msconvert \\
            $raw_file \\
            --outdir . \\
            $peak_picking_filter \\
            $args
        CONVERTER_NAME="msconvert"
        CONVERTER_VERSION=\$(msconvert 2>&1 | grep -oP 'ProteoWizard release: \\K[0-9.]+' || true)
    fi

    # Rename output to include prefix (sample name)
    # msconvert uses the input filename, so we need to rename
    original_output=\$(ls *.mzML | head -1)
    if [ "\$original_output" != "${output_name}" ]; then
        mv "\$original_output" "${output_name}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \${CONVERTER_NAME}: \${CONVERTER_VERSION:-unknown}
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
