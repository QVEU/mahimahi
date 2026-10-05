#!/bin/bash
# End-to-end tests for the mahimahi Snakemake workflow.
#
#   bash tests/workflow/run_workflow_tests.sh
#
# Needs snakemake, pysam, pandas, minimap2 and samtools on PATH
# (workflow/envs/mahimahi.yaml pins them).
#
# Generates synthetic reads whose per-cell, per-template, per-strand UMI
# counts are known exactly, runs the workflow over all three supported input
# types, and checks the output against that truth. Also confirms that an
# inverted strand convention is caught rather than silently producing
# upside-down ratios.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1

REGENERATE=0
REQUIRE_DEPS=0
for arg in "$@"; do
    case "$arg" in
        --regenerate)   REGENERATE=1 ;;
        --require-deps) REQUIRE_DEPS=1 ;;
        *) echo "unknown option: $arg" >&2
           echo "usage: $0 [--regenerate] [--require-deps]" >&2; exit 2 ;;
    esac
done

# ---------------------------------------------------------------------------
# Dependency check
#
# Reports exactly what is missing, and which interpreter was inspected. A bare
# "missing" is not useful on a cluster, where the usual cause is not an absent
# package but the wrong python being found -- see report_python_environment.
# ---------------------------------------------------------------------------
report_python_environment() {
    echo "  python3:  $(command -v python3 || echo '<not found>')" >&2
    if command -v python3 >/dev/null; then
        echo "  version:  $(python3 --version 2>&1)" >&2
    fi
    if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
        echo "  conda env active: ${CONDA_DEFAULT_ENV}" >&2
        if [[ -n "${LOADEDMODULES:-}" ]]; then
            echo >&2
            echo "  NOTE: a conda environment is active AND Lmod modules are loaded." >&2
            echo "  An active conda env puts its own python3 first on PATH, so" >&2
            echo "  'module load py-pysam' has no effect on which interpreter runs." >&2
            echo "  Either 'conda deactivate' before loading modules, or install the" >&2
            echo "  dependencies into a conda env and skip the modules entirely." >&2
        fi
    fi
}

suggest_conda_env() {
    cat >&2 <<'HINT'

  Recommended: one conda environment with all of it, from the file in this repo.

      conda env create -f workflow/envs/mahimahi.yaml
      conda activate mahimahi
      bash tests/workflow/run_workflow_tests.sh

  Or let Snakemake manage it per rule:

      snakemake --cores 8 --software-deployment-method conda

  On an Lmod cluster, note that py-pysam and py-pandas may be built against
  DIFFERENT python versions, in which case loading one unloads the other's
  interpreter and they cannot both be active. Check with:

      module load py-pysam py-pandas && python3 -c 'import pysam, pandas'

  If that fails, the conda route above avoids the conflict.
HINT
}

missing_tools=()
for tool in snakemake minimap2 samtools python3; do
    command -v "$tool" >/dev/null || missing_tools+=("$tool")
done

missing_modules=()
if command -v python3 >/dev/null; then
    for mod in pysam pandas; do
        python3 -c "import ${mod}" 2>/dev/null || missing_modules+=("${mod}")
    done
else
    missing_modules=(pysam pandas)
fi

if [[ ${#missing_tools[@]} -gt 0 || ${#missing_modules[@]} -gt 0 ]]; then
    if [[ "${REQUIRE_DEPS}" -eq 1 ]]; then
        echo "FAIL: dependencies missing and --require-deps was given." >&2
    else
        echo "SKIP: cannot run the workflow tests." >&2
    fi
    [[ ${#missing_tools[@]} -gt 0 ]] &&         echo "  not on PATH:        ${missing_tools[*]}" >&2
    [[ ${#missing_modules[@]} -gt 0 ]] &&         echo "  not importable:     ${missing_modules[*]}" >&2
    echo >&2
    report_python_environment
    # Show the real import error for the first missing module; "missing" is
    # often actually a broken build or an ABI mismatch.
    if [[ ${#missing_modules[@]} -gt 0 ]] && command -v python3 >/dev/null; then
        echo >&2
        echo "  import error for '${missing_modules[0]}':" >&2
        python3 -c "import ${missing_modules[0]}" 2>&1 | sed 's/^/    /' >&2
    fi
    suggest_conda_env
    # A silently skipped suite is worse than a red one in CI, where nobody is
    # watching the output.
    [[ "${REQUIRE_DEPS}" -eq 1 ]] && exit 1
    exit 0
fi

echo "dependencies OK:"
echo "  python3    $(python3 --version 2>&1 | awk '{print $2}')  ($(command -v python3))"
echo "  pysam      $(python3 -c 'import pysam; print(pysam.__version__)')"
echo "  pandas     $(python3 -c 'import pandas; print(pandas.__version__)')"
echo "  snakemake  $(snakemake --version 2>&1 | tail -1)"
echo "  minimap2   $(minimap2 --version 2>&1 | head -1)"
echo "  samtools   $(samtools --version 2>&1 | head -1 | awk '{print $2}')"

CFG=tests/workflow/config_test.yaml
failures=0
step() { echo; echo "### $*"; }

step "1. Example data"
# The fixtures are committed, so a fresh clone can run this immediately. Verify
# they still match the generator rather than trusting that they do; regenerate
# only when asked.
if [[ "${REGENERATE}" -eq 1 ]]; then
    python3 tests/workflow/make_fixtures.py || { echo "FAILED: generation"; exit 1; }
else
    python3 tests/workflow/make_fixtures.py --check \
        || { echo "FAILED: committed fixtures do not match the generator"; failures=$((failures+1)); }
fi

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
MAHIMAHI_REQUIRE_DEPS="${REQUIRE_DEPS}" python3 tests/workflow/verify_against_truth.py || failures=$((failures+1))

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
