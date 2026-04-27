process SQANTI_PROTEIN {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container 'docker://docker.io/anaconesalab/sqanti3:v6.0.1'

    input:
    tuple val(meta), path(cds_gtf)
    path reference_gtf
    path sqanti_protein_script

    output:
    tuple val(meta), path("*.predicted_proteome.best_ORF_SQANTI_classification.tsv"), emit: protein_classification
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    source /conda/miniconda3/etc/profile.d/conda.sh
    conda activate sqanti3
    export SQANTI_PATH=\$(dirname \$(which sqanti3_qc.py))

    # Add src/utilities and utilities to PYTHONPATH for cupcake and other imports
    export PYTHONPATH=\${SQANTI_PATH}/src/utilities:\${SQANTI_PATH}/utilities:\${SQANTI_PATH}:\${PYTHONPATH:-}

    # Copy the script locally and patch it to use the platform-specific gtfToGenePred binary
    # v5.5.4 has gtfToGenePred-linux-x86_64 instead of gtfToGenePred
    cp $sqanti_protein_script ./sqanti3_protein_input_full_gtf_patched.py
    sed -i 's|GTF2GENEPRED_PROG = os.path.join(sqanti_path, "src", "utilities", "gtfToGenePred")|GTF2GENEPRED_PROG = os.path.join(sqanti_path, "src", "utilities", "gtfToGenePred-linux-x86_64")|g' ./sqanti3_protein_input_full_gtf_patched.py

    python ./sqanti3_protein_input_full_gtf_patched.py \\
        $cds_gtf \\
        $reference_gtf \\
        -d . \\
        -p $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.5.4
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.predicted_proteome.best_ORF_SQANTI_classification.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
