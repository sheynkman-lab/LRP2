process SQANTI_PROTEIN {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/anaconesalab/sqanti3:v5.5.4'

    input:
    tuple val(meta), path(cds_gtf)
    path reference_gtf
    path sqanti_protein_script

    output:
    tuple val(meta), path("*.predicted_proteome.best_ORF_SQANTI_classification.tsv"), emit: protein_classification
    tuple val(meta), path("*_S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    exec > >(tee ${prefix}_S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt) 2>&1

    # The anaconesalab/sqanti3 image keeps SQANTI3 in a conda env of its own. Under -profile conda,
    # Nextflow has already activated this module's environment.yml and that path doesn't exist.
    if [ -f /conda/miniconda3/etc/profile.d/conda.sh ]; then
        source /conda/miniconda3/etc/profile.d/conda.sh
        conda activate sqanti3
    fi

    # SQANTI3's install folder. readlink -f follows symlinks: bioconda puts a symlink in <env>/bin
    # that points to <env>/share/sqanti3-<version>-<build>/, where src/utilities lives.
    SQANTI_QC_SCRIPT=\$(command -v sqanti3_qc.py || true)
    if [ -z "\$SQANTI_QC_SCRIPT" ]; then
        echo "ERROR: sqanti3_qc.py not found on PATH; SQANTI3 is missing from this environment" >&2
        exit 1
    fi
    export SQANTI_PATH=\$(dirname "\$(readlink -f "\$SQANTI_QC_SCRIPT")")

    # Add src/utilities and utilities to PYTHONPATH for cupcake and other imports
    export PYTHONPATH=\${SQANTI_PATH}/src/utilities:\${SQANTI_PATH}/utilities:\${SQANTI_PATH}:\${PYTHONPATH:-}

    # The script expects src/utilities/gtfToGenePred. SQANTI3 v5.5.4 ships it as
    # gtfToGenePred-linux-x86_64, and other installs may only have it on PATH (ucsc-gtftogenepred).
    GTF2GENEPRED=""
    for candidate in \\
        "\${SQANTI_PATH}/src/utilities/gtfToGenePred" \\
        "\${SQANTI_PATH}/src/utilities/gtfToGenePred-linux-x86_64" \\
        "\$(command -v gtfToGenePred || true)"; do
        if [ -n "\$candidate" ] && [ -x "\$candidate" ]; then
            GTF2GENEPRED="\$candidate"
            break
        fi
    done
    if [ -z "\$GTF2GENEPRED" ]; then
        echo "ERROR: gtfToGenePred not found in \${SQANTI_PATH}/src/utilities or on PATH" >&2
        exit 1
    fi
    echo "Using gtfToGenePred: \$GTF2GENEPRED"

    # Copy the script locally and point it at the gtfToGenePred found above
    cp $sqanti_protein_script ./sqanti3_protein_input_full_gtf_patched.py
    if [ "\$GTF2GENEPRED" != "\${SQANTI_PATH}/src/utilities/gtfToGenePred" ]; then
        sed -i 's|GTF2GENEPRED_PROG = os.path.join(sqanti_path, "src", "utilities", "gtfToGenePred")|GTF2GENEPRED_PROG = "'"\$GTF2GENEPRED"'"|' \\
            ./sqanti3_protein_input_full_gtf_patched.py
        if ! grep -qF "GTF2GENEPRED_PROG = \\"\$GTF2GENEPRED\\"" ./sqanti3_protein_input_full_gtf_patched.py; then
            echo "ERROR: could not point $sqanti_protein_script at \$GTF2GENEPRED (GTF2GENEPRED_PROG line not found)" >&2
            exit 1
        fi
    fi

    python ./sqanti3_protein_input_full_gtf_patched.py \\
        $cds_gtf \\
        $reference_gtf \\
        -d . \\
        -p $prefix \\
        $args

    # Keep this in step with the container tag above and the version in environment.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.5.4
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.predicted_proteome.best_ORF_SQANTI_classification.tsv
    touch ${prefix}_S3_PREDICTED_PROTEOME_M3_SQANTI_PROTEIN_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.5.4
    END_VERSIONS
    """
}
