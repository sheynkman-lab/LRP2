#!/usr/bin/env bash
#
# Re-create the expected-output snapshots from known-good pipeline results (e.g. HPC runs).
# Run from the root of the LRP2 repository, then commit tests/output_snapshots/*.json.
#
# Usage:
#   tests/output_snapshots/make_snapshots.sh /path/to/LRP2_test_data_expected_outputs_261107
#
# Optional: PIPELINE_DIR=/path/to/the/LRP2/checkout/used/for/those/runs
#   Paths inside that folder (e.g. to bin/ scripts) are written as <PIPELINE_DIR>, the same token
#   the CI comparison uses for its own checkout, so they don't count as differences.
#
set -euo pipefail

EXPECTED_ROOT="${1:?usage: $0 /path/to/expected_outputs_dir}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${HERE}/../../.github/scripts/output_snapshot.py"

# snapshot name (used by integration_tests.yml)  <-  results folder inside EXPECTED_ROOT
SNAPSHOTS=(
    "custom_v42_chr22_test_rna  test_custom_v42_chr22_testrna"
    "test_rna_GRCh38.p14.v50    test_results_test_rna_GRCh38.p14.v50"
    "test_rna_GRCh38.p13.v40    test_results_test_rna_GRCh38.p13.v40"
    "test_dda_GRCh38.p14.v50    test_results_test_dda_GRCh38.p14.v50"
    "test_dia_GRCh38.p14.v50    test_results_test_dia_GRCh38.p14.v50"
    "s2_entry                   test_results_s2_entry"
    "s4_entry                   test_results_s4_entry"
    "s5_entry                   test_results_s5_entry"
)

replace=()
if [[ -n "${PIPELINE_DIR:-}" ]]; then
    replace=(--replace "${PIPELINE_DIR%/}=<PIPELINE_DIR>")
fi

for entry in "${SNAPSHOTS[@]}"; do
    read -r name folder <<< "$entry"
    if [[ ! -d "${EXPECTED_ROOT}/${folder}" ]]; then
        echo "SKIP  ${name}: ${EXPECTED_ROOT}/${folder} not found" >&2
        continue
    fi
    python3 "$SCRIPT" snapshot "${EXPECTED_ROOT}/${folder}" -o "${HERE}/${name}.json" "${replace[@]}"
done
