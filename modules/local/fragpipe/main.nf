/*
 * FRAGPIPE: Comprehensive proteomics analysis with FragPipe
 *
 * This module handles the complete FragPipe workflow:
 * 1. Generate manifest file from mzML files and metadata
 * 2. Prepare protein database (add decoys with philosopher if needed)
 * 3. Download or use custom workflow file (DDA: LFQ-MBR, DIA: DIA_SpecLib_Quant)
 * 4. Execute FragPipe in headless mode
 */

process FRAGPIPE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "https://cloud.sylabs.io/library/jtllab/default/fragpipe:latest"

    input:
    tuple val(meta), path(mzml_files), path(protein_fasta)
    path msfragger_jar
    path ionquant_jar
    path diatracer_jar
    path custom_workflow

    output:
    tuple val(meta), path("*.tsv"), emit: psm_table, optional: true
    tuple val(meta), path("results/**"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def data_type = meta.mass_spec_type ?: 'DDA'
    def replicate = meta.replicate ?: '1'
    def experiment = meta.id
    def decoy_tag = params.fragpipe_decoy_tag ?: 'rev_'
    def mzml_list = mzml_files instanceof List ? mzml_files : [mzml_files]
    def use_custom_workflow = custom_workflow.name != 'NO_FILE'
    def workflow_name = data_type == 'DIA' ? 'DIA_SpecLib_Quant' : 'LFQ-MBR'
    """
    echo "=========================================================================="
    echo "                      FRAGPIPE ANALYSIS: ${meta.id}"
    echo "=========================================================================="
    echo "Data type: ${data_type}"
    echo "Experiment: ${experiment}"
    echo "Replicate: ${replicate}"
    echo "Number of files: ${mzml_list.size()}"
    echo "Database: ${protein_fasta}"
    echo "Decoy tag: ${decoy_tag}"
    echo "=========================================================================="
    echo ""

    #
    # STEP 1: CREATE MANIFEST FILE
    #
    echo "[1/4] Creating FragPipe manifest file..."
    rm -f manifest.fp-manifest

    ${mzml_list.collect { file ->
        "echo \"\$(realpath ${file})\t${experiment}\t${replicate}\t${data_type}\" >> manifest.fp-manifest"
    }.join('\n    ')}

    echo "Manifest contents:"
    cat manifest.fp-manifest
    echo ""

    #
    # STEP 2: PREPARE DATABASE (add decoys if needed)
    #
    echo "[2/4] Preparing protein database..."

    # Check if decoys already exist in the FASTA file
    if grep -q "^>${decoy_tag}" ${protein_fasta}; then
        echo "Decoys found in FASTA file (prefix: ${decoy_tag})"
        cp ${protein_fasta} database_with_decoys.fasta
    else
        echo "No decoys found - generating decoys with philosopher..."
        PHILOSOPHER=\$(find /fragpipe_bin -name "philosopher*" -type f -executable 2>/dev/null | head -1)
        if [ -z "\$PHILOSOPHER" ]; then
            echo "ERROR: Philosopher not found in FragPipe installation"
            exit 1
        fi

        echo "Using philosopher: \$PHILOSOPHER"
        cp ${protein_fasta} reference_proteome.fasta

        \$PHILOSOPHER workspace --init
        \$PHILOSOPHER database --custom reference_proteome.fasta --prefix ${decoy_tag}
        \$PHILOSOPHER workspace --clean

        mv *.fas database_with_decoys.fasta

        echo "Decoys added successfully!"
    fi

    # Report database statistics
    TOTAL_SEQS=\$(grep -c "^>" database_with_decoys.fasta)
    DECOY_SEQS=\$(grep -c "^>${decoy_tag}" database_with_decoys.fasta 2>/dev/null || echo 0)
    DECOY_SEQS=\$(echo "\$DECOY_SEQS" | tr -d '[:space:]')
    TARGET_SEQS=\$((TOTAL_SEQS - DECOY_SEQS))

    echo ""
    echo "Database statistics:"
    echo "  Total sequences: \$TOTAL_SEQS"
    echo "  Target sequences: \$TARGET_SEQS"
    echo "  Decoy sequences: \$DECOY_SEQS"
    echo ""

    #
    # STEP 3: GET WORKFLOW FILE
    #
    echo "[3/4] Preparing workflow file..."

    if [ "${use_custom_workflow}" = "true" ]; then
        echo "Using custom workflow: ${custom_workflow}"
        cp ${custom_workflow} workflow.workflow
    else
        echo "Downloading default ${data_type} workflow (${workflow_name})..."

        # Determine workflow URL
        if [ "${data_type}" = "DIA" ]; then
            WORKFLOW_URL="https://raw.githubusercontent.com/Nesvilab/FragPipe/develop/workflows/DIA_SpecLib_Quant.workflow"
        else
            WORKFLOW_URL="https://raw.githubusercontent.com/Nesvilab/FragPipe/develop/workflows/LFQ-MBR.workflow"
        fi

        # Download workflow (try multiple methods)I
        if command -v wget &> /dev/null; then
            wget -q -O workflow.workflow "\$WORKFLOW_URL"
        elif command -v curl &> /dev/null; then
            curl -sL -o workflow.workflow "\$WORKFLOW_URL"
        elif command -v python3 &> /dev/null; then
            python3 -c "import urllib.request; urllib.request.urlretrieve('\$WORKFLOW_URL', 'workflow.workflow')"
        else
            echo "ERROR: No download tool found (wget, curl, or python)"
            exit 1
        fi

        if [ ! -f workflow.workflow ] || [ ! -s workflow.workflow ]; then
            echo "ERROR: Failed to download workflow file"
            exit 1
        fi

        echo "Workflow downloaded: ${workflow_name}.workflow"
    fi

    # Update workflow with database path and decoy tag
    WORK_DIR="\$(pwd)"
    echo "database.db-path=\$WORK_DIR/database_with_decoys.fasta" >> workflow.workflow
    echo "database.decoy-tag=${decoy_tag}" >> workflow.workflow

    echo ""

    #
    # STEP 4: RUN FRAGPIPE
    #
    echo "[4/4] Running FragPipe..."
    echo ""

    # Create writable directories needed for FragPipe
    # FragPipe expects: config_tools/msfragger/*.jar, config_tools/ionquant/*.jar, config_tools/diatracer/*.jar
    mkdir -p results fragpipe_tools/msfragger fragpipe_tools/ionquant fragpipe_tools/diatracer
    cp ${msfragger_jar} fragpipe_tools/msfragger/
    cp ${ionquant_jar} fragpipe_tools/ionquant/
    cp ${diatracer_jar} fragpipe_tools/diatracer/
    chmod +x fragpipe_tools/msfragger/*.jar fragpipe_tools/ionquant/*.jar fragpipe_tools/diatracer/*.jar

    # Set working directory for absolute paths in manifest/workflow
    WORK_DIR="\$(pwd)"

    echo "FragPipe tools:"
    echo "  MSFragger: \$(ls -lh fragpipe_tools/msfragger/*.jar | awk '{print \$9, \$5}')"
    echo "  IonQuant: \$(ls -lh fragpipe_tools/ionquant/*.jar | awk '{print \$9, \$5}')"
    echo "  DiaTracer: \$(ls -lh fragpipe_tools/diatracer/*.jar | awk '{print \$9, \$5}')"
    echo ""

    # Find FragPipe executable
    if command -v fragpipe &> /dev/null; then
        FRAGPIPE_CMD="fragpipe"
    elif [ -f "/fragpipe_bin/fragpipe-${params.fragpipe_version}/fragpipe-${params.fragpipe_version}/bin/fragpipe" ]; then
        FRAGPIPE_CMD="/fragpipe_bin/fragpipe-${params.fragpipe_version}/fragpipe-${params.fragpipe_version}/bin/fragpipe"
    elif [ -d "/fragpipe_bin" ]; then
        FRAGPIPE_CMD=\$(find /fragpipe_bin -name "fragpipe" -type f -executable 2>/dev/null | head -1)
        if [ -z "\$FRAGPIPE_CMD" ]; then
            echo "ERROR: FragPipe executable not found in /fragpipe_bin"
            exit 1
        fi
    else
        echo "ERROR: FragPipe executable not found"
        exit 1
    fi

    echo "FragPipe executable: \$FRAGPIPE_CMD"
    echo "RAM: ${task.memory.toGiga()} GB"
    echo "Threads: ${task.cpus}"
    echo ""
    echo "Starting FragPipe headless execution..."
    echo "----------------------------------------------------------------------"

    \$FRAGPIPE_CMD --headless \\
        --workflow \$WORK_DIR/workflow.workflow \\
        --manifest \$WORK_DIR/manifest.fp-manifest \\
        --workdir \$WORK_DIR/results \\
        --config-tools-folder \$WORK_DIR/fragpipe_tools \\
        --ram ${task.memory.toGiga()} \\
        --threads ${task.cpus} \\
        $args

    EXIT_CODE=\$?

    echo "----------------------------------------------------------------------"

    if [ \$EXIT_CODE -ne 0 ]; then
        echo ""
        echo "ERROR: FragPipe failed with exit code \$EXIT_CODE"
        echo ""
        echo "Diagnostic information:"
        echo "========================================================================"
        echo "Results directory contents:"
        ls -lah results/ || true
        echo ""
        echo "========================================================================"
        echo "FragPipe log (last 100 lines):"
        find results -name "log*.txt" -exec tail -100 {} \\; 2>/dev/null || echo "No log file found"
        echo "========================================================================"
        exit 1
    fi

    echo ""
    echo "FragPipe completed successfully!"
    echo ""

    #
    # POST-PROCESSING: Extract PSM table
    #
    echo "Searching for output files..."

    # Try different PSM table names (depends on workflow type)
    if find results -name "psm.tsv" -exec cp {} ${prefix}.tsv \\; 2>/dev/null; then
        echo "PSM table found: psm.tsv"
    elif find results -name "*_psm.tsv" -exec cp {} ${prefix}.tsv \\; 2>/dev/null; then
        echo "PSM table found: *_psm.tsv"
    elif find results -name "combined_peptide.tsv" -exec cp {} ${prefix}.tsv \\; 2>/dev/null; then
        echo "Output table found: combined_peptide.tsv"
    elif find results -name "combined_protein.tsv" -exec cp {} ${prefix}.tsv \\; 2>/dev/null; then
        echo "Output table found: combined_protein.tsv"
    else
        echo "Note: No standard PSM table found."
        echo "All results are available in the results directory."
    fi

    echo ""
    echo "Results directory structure:"
    find results -maxdepth 2 -type f | head -20
    echo ""
    echo "=========================================================================="
    echo "FragPipe analysis completed for: ${meta.id}"
    echo "=========================================================================="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragpipe: \$(\$FRAGPIPE_CMD --version 2>&1 | grep -oP 'FragPipe \\K[0-9.]+' || echo "${params.fragpipe_version}")
        msfragger: \$(basename ${msfragger_jar} | grep -oP 'MSFragger-\\K[0-9.]+' || echo "unknown")
        ionquant: \$(basename ${ionquant_jar} | grep -oP 'IonQuant-\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results
    touch ${prefix}.tsv
    touch results/psm.tsv
    touch results/combined_protein.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragpipe: ${params.fragpipe_version}
        msfragger: 4.1
        ionquant: 1.10.12
    END_VERSIONS
    """
}
