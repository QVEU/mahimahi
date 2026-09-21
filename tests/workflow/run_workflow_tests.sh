#!/bin/bash
# End-to-end tests for the SCISSORS Snakemake workflow.
#
#   bash tests/workflow/run_workflow_tests.sh
#
# Needs snakemake, pysam, pandas, minimap2 and samtools on PATH
# (workflow/envs/scissors.yaml pins them).
#
# Generates synthetic reads whose per-cell, per-template, per-strand UMI
# counts are known exactly, runs the workflow over all three supported input
# types, and checks the output against that truth. Also confirms that an
# inverted strand convention is caught rather than silently producing
# upside-down ratios.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1

for tool in snakemake minimap2 samtools python3; do
    command -v "$tool" >/dev/null || { echo "SKIP: $tool not on PATH" >&2; exit 0; }
done
python3 -c 'import pysam, pandas' 2>/dev/null || { echo "SKIP: pysam/pandas missing" >&2; exit 0; }

CFG=tests/workflow/config_test.yaml
failures=0
step() { echo; echo "### $*"; }

step "1. Generate fixtures with known ground truth"
python3 tests/workflow/make_fixtures.py || { echo "FAILED"; exit 1; }

step "2. Lint the workflow"
snakemake --configfile "$CFG" --cores 1 --lint 2>&1 | tail -20

step "3. Dry run"
snakemake --configfile "$CFG" --cores 4 -n >/dev/null 2>&1 \
    && echo "dry run OK" || { echo "FAILED: dry run"; failures=$((failures+1)); }

step "4. Run the workflow over FASTQ+DRAGEN, raw 10x paired FASTQ, and tagged BAM"
rm -rf tests/workflow/results
if snakemake --configfile "$CFG" --cores 4 2>&1 | tail -3; then
    echo "workflow completed"
else
    echo "FAILED: workflow run"; failures=$((failures+1))
fi

step "5. Verify counts against ground truth"
python3 tests/workflow/verify_against_truth.py || failures=$((failures+1))

step "6. Confirm an inverted strand convention is rejected"
sed -e 's/convention: "reverse_is_positive"/convention: "forward_is_positive"/' \
    -e 's#results_dir: "tests/workflow/results"#results_dir: "tests/workflow/results_inv"#' \
    "$CFG" > tests/workflow/config_inverted.yaml
rm -rf tests/workflow/results_inv
if snakemake --configfile tests/workflow/config_inverted.yaml --cores 4 >/dev/null 2>&1; then
    echo "FAILED: inverted convention was accepted; the sense control should have caught it"
    failures=$((failures+1))
else
    if grep -q "convention is probably inverted" \
         tests/workflow/results_inv/*/logs/strand_qc.log 2>/dev/null; then
        echo "inverted convention correctly rejected by the sense control"
    else
        echo "FAILED: run failed, but not via the strand_qc sense control"
        failures=$((failures+1))
    fi
fi
rm -rf tests/workflow/results_inv tests/workflow/config_inverted.yaml

step "7. Idempotence: a second run should do nothing"
out=$(snakemake --configfile "$CFG" --cores 4 2>&1)
if echo "$out" | grep -qE 'Nothing to be done|nothing to be done'; then
    echo "second run is a no-op, as expected"
else
    echo "FAILED: re-run was not a no-op"; failures=$((failures+1))
fi

echo
echo "################################################################"
[[ "$failures" -eq 0 ]] && echo "# Workflow tests: all passed." \
                        || echo "# Workflow tests: $failures FAILED."
echo "################################################################"
[[ "$failures" -eq 0 ]]
