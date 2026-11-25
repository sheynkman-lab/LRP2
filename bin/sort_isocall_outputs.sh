#!/usr/bin/env bash
#
# Rewrite isocall call's outputs in a fixed order, so the same input always gives
# byte-identical files (isocall call writes isoforms in whatever order its threads finish).
#
#   sort_isocall_outputs.sh PREFIX.isocall.isoforms.gtf.gz PREFIX.isocall.count_matrix.csv
#
# GTF: comment lines first, then transcripts ordered by chromosome, transcript start and
#      transcript_id; each transcript's own lines stay together (transcript line, then exons
#      by position, then any other features).
# Count matrix: header line first, then rows sorted by transcript ID.
# Only the order changes: the files keep exactly the same lines.
set -euo pipefail
export LC_ALL=C

gtf="${1:?usage: $0 ISOFORMS.gtf.gz COUNT_MATRIX.csv}"
matrix="${2:?usage: $0 ISOFORMS.gtf.gz COUNT_MATRIX.csv}"
tmp="$(mktemp -d "${TMPDIR:-.}/sort_isocall.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

# --- GTF ---------------------------------------------------------------------------------
gzip -cd "$gtf" > "$tmp/in.gtf"
n_in=$(wc -l < "$tmp/in.gtf")
{
    grep '^#' "$tmp/in.gtf" || true
    # Pass 1 records each transcript's start; pass 2 prefixes every line with sort keys:
    #   chrom, transcript start, transcript_id, feature rank, start, end, original line
    awk -F '\t' -v OFS='\t' '
        function tid(attr) {
            if (match(attr, /transcript_id "[^"]*"/)) return substr(attr, RSTART + 15, RLENGTH - 16)
            return ""
        }
        FNR == NR {
            if ($0 !~ /^#/ && NF >= 9) {
                t = tid($9)
                if (t != "" && (!(t in start) || $4 + 0 < start[t])) start[t] = $4 + 0
            }
            next
        }
        /^#/ { next }
        {
            t = tid($9)
            s = (t in start) ? start[t] : $4 + 0
            rank = ($3 == "gene") ? 0 : ($3 == "transcript") ? 1 : ($3 == "exon") ? 2 : 3
            print $1, s, t, rank, $4 + 0, $5 + 0, $0
        }' "$tmp/in.gtf" "$tmp/in.gtf" \
    | sort -t "$(printf '\t')" -k1,1 -k2,2n -k3,3 -k4,4n -k5,5n -k6,6n -k7 \
    | cut -f 7-
} > "$tmp/out.gtf"
n_out=$(wc -l < "$tmp/out.gtf")
if [ "$n_in" -ne "$n_out" ]; then
    echo "sort_isocall_outputs.sh: GTF line count changed ($n_in -> $n_out); leaving $gtf unchanged" >&2
    exit 1
fi
gzip -n -c "$tmp/out.gtf" > "$tmp/out.gtf.gz"
mv "$tmp/out.gtf.gz" "$gtf"

# --- Count matrix ------------------------------------------------------------------------
{
    head -n 1 "$matrix"
    tail -n +2 "$matrix" | sort -t ',' -k1,1 -k2
} > "$tmp/matrix"
if [ "$(wc -l < "$matrix")" -ne "$(wc -l < "$tmp/matrix")" ]; then
    echo "sort_isocall_outputs.sh: count matrix line count changed; leaving $matrix unchanged" >&2
    exit 1
fi
cat "$tmp/matrix" > "$matrix"
