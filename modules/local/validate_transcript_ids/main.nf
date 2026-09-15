process VALIDATE_TRANSCRIPT_IDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker://ubuntu:22.04"

    input:
    tuple val(meta), path(gtf), path(counts)

    output:
    tuple val(meta), path(gtf), path(counts), emit: validated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export TMPDIR=\$(pwd)/tmp
    mkdir -p \$TMPDIR

    if [[ "$gtf" == *.gz ]]; then
        zcat "$gtf" | grep -v "^#" | grep -oP 'transcript_id "\\K[^"]+' | sort -u > gtf_ids.txt
    else
        grep -v "^#" "$gtf" | grep -oP 'transcript_id "\\K[^"]+' | sort -u > gtf_ids.txt
    fi

    tail -n +2 "$counts" | cut -f1 | sort -u > counts_ids.txt
    GTF_COUNT=\$(wc -l < gtf_ids.txt)
    COUNTS_COUNT=\$(wc -l < counts_ids.txt)

    comm -12 gtf_ids.txt counts_ids.txt > matching_ids.txt
    MATCHING_COUNT=\$(wc -l < matching_ids.txt)

    if [ \$COUNTS_COUNT -gt 0 ]; then
        OVERLAP_PCT=\$(awk "BEGIN {printf \\"%.1f\\", (\$MATCHING_COUNT / \$COUNTS_COUNT) * 100}")
    else
        OVERLAP_PCT=0
    fi

    echo "Transcript ID Validation Results:"
    echo "  GTF unique transcript IDs: \$GTF_COUNT"
    echo "  Counts matrix unique transcript IDs: \$COUNTS_COUNT"
    echo "  Matching transcript IDs: \$MATCHING_COUNT"
    echo "  Overlap: \${OVERLAP_PCT}%"

    # Fail and exit if no overlap: users should correct their own data inputs.
    if [ \$MATCHING_COUNT -eq 0 ]; then
        echo ""
        echo "ERROR: No matching transcript IDs found between GTF and count matrix!"
        echo ""
        echo "GTF IDs (first 10):"
        head -10 gtf_ids.txt
        echo ""
        echo "Counts IDs (first 10):"
        head -10 counts_ids.txt
        echo ""
        echo "Please ensure that your inputted GTF and counts matrix files are from the same isoform quantification run."
        exit 1
    fi

    comm -23 gtf_ids.txt counts_ids.txt > only_in_gtf.txt
    comm -13 gtf_ids.txt counts_ids.txt > only_in_counts.txt

    ONLY_GTF_COUNT=\$(wc -l < only_in_gtf.txt)
    ONLY_COUNTS_COUNT=\$(wc -l < only_in_counts.txt)

    # Fail if there are ANY mismatched transcripts
    if [ \$ONLY_GTF_COUNT -gt 0 ] || [ \$ONLY_COUNTS_COUNT -gt 0 ]; then
        echo ""
        echo "ERROR: Transcript ID mismatch detected between GTF and counts matrix!"
        echo ""
        echo "Summary:"
        echo "  GTF unique transcript IDs: \$GTF_COUNT"
        echo "  Counts matrix unique transcript IDs: \$COUNTS_COUNT"
        echo "  Matching transcript IDs: \$MATCHING_COUNT"
        echo "  Transcripts ONLY in GTF (not in counts): \$ONLY_GTF_COUNT"
        echo "  Transcripts ONLY in counts (not in GTF): \$ONLY_COUNTS_COUNT"
        echo ""

        if [ \$ONLY_GTF_COUNT -gt 0 ]; then
            echo "Transcripts in GTF but NOT in counts (first 20):"
            head -20 only_in_gtf.txt
            echo ""
        fi

        if [ \$ONLY_COUNTS_COUNT -gt 0 ]; then
            echo "Transcripts in counts but NOT in GTF (first 20):"
            head -20 only_in_counts.txt
            echo ""
        fi

        echo "  WARNING: Your inputted GTF and counts files must have exactly the same transcript IDs."
        echo "  Please filter your input files to only include matching transcripts."
        exit 1
    fi

    echo ""
    echo "All transcript IDs match perfectly between GTF and counts!"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -1 | awk '{print \$4}')
    END_VERSIONS
    """

    stub:
    """
    touch gtf_ids.txt counts_ids.txt matching_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -1 | awk '{print \$4}')
    END_VERSIONS
    """
}
