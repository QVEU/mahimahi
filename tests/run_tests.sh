#!/bin/bash
# Offline regression tests for the SCISSORS pipeline.
#
#   bash tests/run_tests.sh
#
# These exercise the correctness fixes -- barcode-keyed doublet joining,
# deterministic mixture component ordering, and Slurm-derived Cell Ranger
# memory -- without needing Cell Ranger, real Seurat, or cluster data.
#
# The R tests need Seurat and mixtools to be loadable. Where the real packages
# are unavailable (no CRAN access), tests/stubs/ holds minimal stand-ins that
# reproduce only the behaviour under test: AddMetaData's named-vs-unnamed
# vector contract, and normalmixEM returning its components in arbitrary
# order. Install them into a throwaway library with --stubs.
#
# Run against the real packages whenever you can; the stubs exist so these
# tests are runnable on a machine with no package access at all.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1

USE_STUBS=0
REQUIRE_DEPS=0
for arg in "$@"; do
    case "$arg" in
        --stubs)        USE_STUBS=1 ;;
        --require-deps) REQUIRE_DEPS=1 ;;
        *) echo "unknown option: $arg" >&2
           echo "usage: $0 [--stubs] [--require-deps]" >&2; exit 2 ;;
    esac
done

if [[ "${USE_STUBS}" -eq 1 ]]; then
    STUB_LIB="$(mktemp -d)"
    trap 'rm -rf "${STUB_LIB}"' EXIT
    echo "Installing test stubs into ${STUB_LIB}"
    for pkg in tests/stubs/*/; do
        R CMD INSTALL --no-docs --no-byte-compile --library="${STUB_LIB}" "${pkg}" >/dev/null 2>&1 \
            || { echo "Failed to install stub ${pkg}" >&2; exit 1; }
    done
    export R_LIBS="${STUB_LIB}:${R_LIBS:-}"
    echo "WARNING: running against stubs, not real Seurat/mixtools." >&2
fi

failures=0
run() {
    echo ""
    echo "################################################################"
    echo "# $1"
    echo "################################################################"
    shift
    if "$@"; then echo "--- OK"; else echo "--- FAILED"; failures=$((failures + 1)); fi
}

run "Slurm memory derivation (config.sh)" bash tests/test_config_memory.sh
run "Shell script syntax" bash -c 'for f in config.sh scripts/*.sh; do bash -n "$f" || exit 1; echo "$f OK"; done'
run "Python script syntax" python3 -c 'import ast; ast.parse(open("scripts/04_scrublet.py").read()); print("04_scrublet.py OK")'
run "R script parsing" Rscript -e 'for (f in c("analysis/config.R","analysis/helpers.R","analysis/PV_mutants_integrated.R","analysis/replication.R","analysis/scissors_replication.R")) { invisible(parse(f)); cat(f,"OK\n") }'
run "R Markdown chunk parsing" Rscript tests/check_rmd_parses.R analysis/CVB3_QC.Rmd analysis/EVA71_QC.Rmd
run "Doublet score joining (helpers.R)" Rscript tests/test_doublet_join.R
run "Infected status calling (helpers.R)" Rscript tests/test_infected_status.R
run "Replication slope fitting (replication.R)" Rscript tests/test_replication_fit.R
if [[ "${USE_STUBS}" -eq 1 ]]; then
    echo
    echo "################################################################"
    echo "# Seurat analyses end-to-end"
    echo "################################################################"
    echo "SKIPPED: --stubs masks the real Seurat with tests/stubs/SeuratStub."
    echo "  These analyses must run against the real package. Re-run without"
    echo "  --stubs, or: bash tests/seurat/run_seurat_smoke.sh"
else
    if [[ "${REQUIRE_DEPS}" -eq 1 ]]; then
        run "Seurat analyses end-to-end" bash tests/seurat/run_seurat_smoke.sh --require-deps
    else
        run "Seurat analyses end-to-end" bash tests/seurat/run_seurat_smoke.sh
    fi
fi

echo ""
echo "################################################################"
if [[ "${failures}" -eq 0 ]]; then
    echo "# All test groups passed."
else
    echo "# ${failures} test group(s) FAILED."
fi
echo "################################################################"
[[ "${failures}" -eq 0 ]]
