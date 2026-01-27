process GENERATE_HASHIDS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(corrected_gtf)
    path hashlib_script
    path generate_hashids_script

    output:
    tuple val(meta), path("sqanti_transcript/*_hashids_mapping.txt"), emit: hashids_mapping
    tuple val(meta), path("sqanti_transcript/*_corrected.psl"), emit: corrected_psl
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Create directory structure expected by the R script
    mkdir -p sqanti_transcript

    # Link the corrected GTF file with expected naming
    ln -sf \$(pwd)/$corrected_gtf sqanti_transcript/${prefix}_corrected.gtf

    # Export environment variables expected by the R script
    export OUTPUT_DIR=\$(pwd)
    export OUTPUT_BASE_NAME=$prefix
    export HASHLIB_SCRIPT=\$(pwd)/$hashlib_script

    # Ensure R can find packages in the container
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    # Run the hash ID generation script
    Rscript \$(pwd)/$generate_hashids_script $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p sqanti_transcript
    touch sqanti_transcript/${prefix}_hashids_mapping.txt
    touch sqanti_transcript/${prefix}_corrected.psl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
        python: 3.11.0
    END_VERSIONS
    """
}
