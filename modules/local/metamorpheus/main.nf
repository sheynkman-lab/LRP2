process METAMORPHEUS {
    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/smithchemwisc/metamorpheus:latest"

    // NOTE: MetaMorpheus needs to write to /MetaMorpheus/CustomAminoAcids and /MetaMorpheus/Mods. 
    // Here we set up all required writable directories and bind mounts for this to work when using container environment. 
    beforeScript """
        mkdir -p \${PWD}/mm_writable_custom_aa \${PWD}/mm_writable_mods \${PWD}/mm_writable_settings
        touch \${PWD}/mm_writable_settings/settings.toml

        # Auto-extract MetaMorpheus default files from container on first run
        METAMORPH_DATA_DIR="${projectDir}/assets/metamorpheus_data"
        CONTAINER="\${NXF_SINGULARITY_CACHEDIR:-${projectDir}/work/singularity}/docker.io-smithchemwisc-metamorpheus-latest.img"

        if [ ! -f "\${METAMORPH_DATA_DIR}/Mods/ProteaseMods.txt" ]; then
            LOCK_FILE="\${METAMORPH_DATA_DIR}/.extraction.lock"
            mkdir -p "\${METAMORPH_DATA_DIR}"

            for i in {1..60}; do
                if mkdir "\${LOCK_FILE}" 2>/dev/null; then
                    echo "Extracting MetaMorpheus default files from container (first-time setup)..."
                    if [ -f "\$CONTAINER" ]; then
                        singularity exec "\$CONTAINER" sh -c 'cd /MetaMorpheus && tar czf - Mods' | \
                            tar xzf - -C "\${METAMORPH_DATA_DIR}/" 2>/dev/null || true

                        singularity exec "\$CONTAINER" sh -c 'cd /MetaMorpheus && tar czf - CustomAminoAcids' | \
                            tar xzf - -C "\${METAMORPH_DATA_DIR}/" 2>/dev/null || true

                        echo "MetaMorpheus files extracted to \${METAMORPH_DATA_DIR}"
                    else
                        echo "Warning: Container not found at \$CONTAINER"
                    fi

                    rmdir "\${LOCK_FILE}"
                    break
                else
                    sleep 1
                fi
            done
        fi

        if [ -d "\${METAMORPH_DATA_DIR}/Mods" ]; then
            cp -r "\${METAMORPH_DATA_DIR}/Mods/"* \${PWD}/mm_writable_mods/ 2>/dev/null || true
        fi

        if [ -d "\${METAMORPH_DATA_DIR}/CustomAminoAcids" ]; then
            cp -r "\${METAMORPH_DATA_DIR}/CustomAminoAcids/"* \${PWD}/mm_writable_custom_aa/ 2>/dev/null || true
        fi
    """
    // Bind mount writable host directories to read-only container paths
    containerOptions '--bind ${PWD}/mm_writable_custom_aa:/MetaMorpheus/CustomAminoAcids --bind ${PWD}/mm_writable_mods:/MetaMorpheus/Mods --bind ${PWD}/mm_writable_settings/settings.toml:/MetaMorpheus/settings.toml'

    input:
    tuple val(meta), path(mzml_files)
    path database_fasta
    path config_toml
    path "mm_writable"

    output:
    tuple val(meta), path("*.psmtsv"), emit: psm_table
    tuple val(meta), path("results/**"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Rename mzML files to simple names based on meta.id
    def mzml_array = mzml_files instanceof List ? mzml_files : [mzml_files]
    def symlink_commands = []
    def renamed_files = []
    mzml_array.eachWithIndex { file, idx ->
        def simple_name = mzml_array.size() > 1 ? "${prefix}_${idx + 1}.mzML" : "${prefix}.mzML"
        symlink_commands.add("ln -s ${file} ${simple_name}")
        renamed_files.add(simple_name)
    }
    def mzml_list = renamed_files.join(' ')

    """
    # Rename mzML files to simple names for cleaner MetaMorpheus output
    ${symlink_commands.join('\n    ')}

    dotnet /MetaMorpheus/CMD.dll \\
        -t $config_toml \\
        -s $mzml_list \\
        -d $database_fasta \\
        -o results/ \\
        $args

    find results -name "AllPeptides*.psmtsv" -exec cp {} ${prefix}.psmtsv \\;
    if [ ! -f "${prefix}.psmtsv" ]; then
        find results -name "*.psmtsv" -type f -print -quit -exec cp {} ${prefix}.psmtsv \\;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamorpheus: \$(dotnet /MetaMorpheus/CMD.dll --version 2>&1 | grep -oP 'MetaMorpheus version \\K[0-9.]+' || echo "1.0.5")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results
    touch ${prefix}.psmtsv
    touch results/AllPeptides.psmtsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamorpheus: 1.0.5
    END_VERSIONS
    """
}
