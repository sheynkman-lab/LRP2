process SANITIZE_SAMPLE_NAMES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker://ubuntu:22.04"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*_sanitized.tsv"), emit: sanitized_counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_name = "${prefix}_samplenames_sanitized.tsv"
    """
    echo "Sanitizing sample column names..."

    # Replace underscores with dashes in sample name columns (header row only)
    awk 'BEGIN {FS=OFS="\\t"}
    NR==1 {
        # Header row: keep first column as-is (transcript_id)
        printf \$1
        # Replace underscores with dashes in sample name columns
        for(i=2; i<=NF; i++) {
            gsub(/_/, "-", \$i)
            printf "\\t" \$i
        }
        printf "\\n"
        next
    }
    {print}' "$counts" > "${output_name}"

    echo "Sample name sanitization summary:"
    SAMPLE_BEFORE=\$(head -1 "$counts" | cut -f2)
    SAMPLE_AFTER=\$(head -1 "${output_name}" | cut -f2)
    if [ "\$SAMPLE_BEFORE" != "\$SAMPLE_AFTER" ]; then
        echo "  Sanitized: '\$SAMPLE_BEFORE' -> '\$SAMPLE_AFTER'"
        echo "  Columns sanitized: \$(head -1 "${output_name}" | awk -F'\\t' '{print NF-1}')"
    else
        echo "  All columns already devoid of underscores, proceeding with remaining modules."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -1 | awk '{print \$4}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_name = "${prefix}_samplenames_sanitized.tsv"
    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -1 | awk '{print \$4}')
    END_VERSIONS
    """
}
