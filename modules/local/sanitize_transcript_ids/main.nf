process SANITIZE_TRANSCRIPT_IDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker://ubuntu:22.04"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*_sanitized.*"), emit: sanitized
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Determine output file extension based on input
    def extension = input_file.name.tokenize('.').last()
    def output_name = extension == "gtf" ? "${input_file.baseName}_sanitized.gtf" : "${prefix}_sanitized.tsv"

    """
    # Sanitize transcript IDs: Replace _PAR_Y with -PAR-Y for pseudoautosomal region gene ID naming in earlier GENCODE refs v25 - v43
    if [[ "$input_file" == *.gtf* ]]; then
        echo "Processing GTF file: $input_file"
        sed 's/_PAR_Y/-PAR-Y/g' "$input_file" > "$output_name"
        CHANGES=\$(grep -c "_PAR_Y" "$input_file" || true)
        if [ \$CHANGES -gt 0 ]; then
            echo "  Sanitized \$CHANGES _PAR_Y occurrences to -PAR-Y"
        else
            echo "  No _PAR_Y IDs found"
        fi

    else
        echo "Processing counts/TSV file: $input_file"
        awk 'BEGIN {FS=OFS="\\t"}
        {
            # Replace _PAR_Y with -PAR-Y in first column (transcript IDs)
            gsub(/_PAR_Y/, "-PAR-Y", \$1)
            print
        }' "$input_file" > "$output_name"

        CHANGES=\$(cut -f1 "$input_file" | grep -c "_PAR_Y" || true)
        if [ \$CHANGES -gt 0 ]; then
            echo "  Sanitized \$CHANGES transcript IDs with _PAR_Y to -PAR-Y"
        else
            echo "  No _PAR_Y IDs found in transcript column"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -1 | awk '{print \$4}')
        sed: \$(sed --version 2>&1 | head -n1 | sed 's/.*) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = input_file.name.tokenize('.').last()
    def output_name = extension == "gtf" ? "${input_file.baseName}_sanitized.gtf" : "${prefix}_sanitized.tsv"
    """
    touch $output_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -1 | awk '{print \$4}')
        sed: \$(sed --version 2>&1 | head -n1 | sed 's/.*) //')
    END_VERSIONS
    """
}
