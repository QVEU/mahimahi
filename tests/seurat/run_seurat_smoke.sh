#!/bin/bash
# Smoke test for the Seurat host-response analyses.
#
#   bash tests/seurat/run_seurat_smoke.sh [--require-deps]
#
# Until this existed, analysis/PV_mutants_integrated.R and the two QC notebooks
# had never executed a single line -- tests/run_tests.sh only parse()d them.
# This synthesises a 10x share layout and runs all three for real.
#
# It found three bugs the parse check could not:
#   - Seurat 5 returns a numeric vector from PercentageFeatureSet where v4
#     returned a one-column data.frame, breaking safe_feature_percentage()'s
#     `[, 1]` with "incorrect number of dimensions".
#   - PV_mutants_integrated.R assigned into a zero-row FindAllMarkers result,
#     failing with "replacement has 1 row, data has 0".
#   - Both notebooks grouped a zero-row marker frame by `cluster`, failing with
#     "Column `cluster` is not found".
#
# --require-deps makes a missing dependency a failure instead of a skip, for
# CI, where a silently skipped suite is worse than a red one.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1

REQUIRE_DEPS=0
[[ "${1:-}" == "--require-deps" ]] && REQUIRE_DEPS=1

skip_or_fail() {
    echo "$1" >&2
    if [[ "$REQUIRE_DEPS" -eq 1 ]]; then
        echo "  (--require-deps given, so this is a failure)" >&2
        exit 1
    fi
    exit 0
}

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
command -v Rscript >/dev/null || skip_or_fail "SKIP: Rscript not on PATH."

missing=$(Rscript -e '
need <- c("Seurat", "SeuratObject", "Matrix", "mixtools", "ggplot2",
          "ggrepel", "dplyr", "patchwork", "rmarkdown", "rprojroot")
absent <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
cat(paste(absent, collapse = " "))' 2>/dev/null)

if [[ -n "$missing" ]]; then
    skip_or_fail "SKIP: R packages not installed: $missing
  The analyses need Seurat 5. Install the environment with:
      conda env create -f workflow/envs/scissors.yaml && conda activate scissors"
fi

seurat_major=$(Rscript -e 'cat(as.integer(packageVersion("Seurat")[1,1]))' 2>/dev/null)
if [[ "${seurat_major:-0}" -lt 5 ]]; then
    skip_or_fail "SKIP: Seurat $(Rscript -e 'cat(as.character(packageVersion("Seurat")))') found, but these analyses require >= 5.0.0.
  Seurat 4 is archived on CRAN; see NEWS.md for what the v5 port changed."
fi

command -v pandoc >/dev/null || skip_or_fail "SKIP: pandoc not on PATH (needed to knit the .Rmd notebooks)."

echo "dependencies OK:"
Rscript -e 'for (p in c("Seurat","SeuratObject","Matrix","mixtools","rmarkdown"))
  cat(sprintf("  %-14s %s\n", p, as.character(packageVersion(p))))' 2>/dev/null
echo "  pandoc         $(pandoc --version | head -1 | awk '{print $2}')"

# ---------------------------------------------------------------------------
ROOT=$(mktemp -d)
OUTDIR=$(mktemp -d)
trap 'rm -rf "$ROOT" "$OUTDIR"' EXIT
export SCISSORS_SHARE_ROOT="$ROOT"

REPO=$PWD
failures=0
step() { echo; echo "### $*"; }
check() {  # description, condition-already-evaluated exit status
    if [[ "$2" -eq 0 ]]; then echo "  PASS  $1"
    else echo "  FAIL  $1"; failures=$((failures + 1)); fi
}

step "1. Synthesise a 10x share layout"
if ! Rscript tests/seurat/make_seurat_fixtures.R "$ROOT" 80 2>&1 | tail -8; then
    echo "FAILED: fixture generation"; exit 1
fi

CR="$ROOT/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger"
RES="$ROOT/Projects/CM_kb"

step "2. Run analysis/PV_mutants_integrated.R (3 of 14 samples)"
export SCISSORS_PV_SAMPLE_IDS="Mock_5h_PV,WT_GFP_PV,RFP_C109S_PV"
( cd "$ROOT" && Rscript "$REPO/analysis/PV_mutants_integrated.R" ) \
    > "$OUTDIR/pv.log" 2>&1
pv_status=$?
check "exits 0" "$pv_status"
[[ "$pv_status" -ne 0 ]] && tail -20 "$OUTDIR/pv.log"

for f in mutsFinal.rds scissors.metadata0.5.csv markers.0.5.csv; do
    [[ -s "$CR/$f" ]]; check "wrote $f" $?
done

step "3. Check the integrated results are sane"
export CR_DIR="$CR"
Rscript -e '
cr <- Sys.getenv("CR_DIR")
md <- read.csv(file.path(cr, "scissors.metadata0.5.csv"), row.names = 1)
mk <- read.csv(file.path(cr, "markers.0.5.csv"))
fail <- 0
say <- function(ok, msg) { cat(if (ok) "  PASS  " else "  FAIL  ", msg, "\n", sep="")
                           if (!ok) fail <<- fail + 1 }

say(nrow(md) > 0, sprintf("%d cells survived QC", nrow(md)))
say(all(c("percent.virus","InfectedStatus","InfectedStatus_groups","doublet_scores",
          "S.Score","G2M.Score","seurat_clusters") %in% colnames(md)),
    "metadata carries viral load, infected status, doublet scores and cell cycle")

## The mock sample takes the uninfected branch and must have no infected cells;
## if the mixture ever ran on it, this is where that shows up.
mock <- md[md$orig.ident == "Mock_5h_PV", ]
say(nrow(mock) > 0 && all(mock$InfectedStatus == "Not_Infected"),
    sprintf("mock: 0/%d cells called infected", nrow(mock)))

## The infected samples must have both classes, or the mixture is degenerate.
for (s in c("WT_GFP_PV", "RFP_C109S_PV")) {
  d <- md[md$orig.ident == s, ]
  n_inf <- sum(d$InfectedStatus == "Infected")
  say(n_inf > 0 && n_inf < nrow(d),
      sprintf("%s: %d/%d infected, both classes present", s, n_inf, nrow(d)))
  say(all(c("High","Low") %in% d$InfectedStatus_groups),
      sprintf("%s: High and Low replication groups both populated", s))
}

say(nrow(mk) > 0 && "cluster" %in% colnames(mk),
    sprintf("FindAllMarkers returned %d rows across %d clusters",
            nrow(mk), length(unique(mk$cluster))))
quit(status = if (fail > 0) 1 else 0)' 2>/dev/null
check "integrated results sane" $?

step "4. Knit the two QC notebooks"
for nb in CVB3_QC EVA71_QC; do
    Rscript -e "rmarkdown::render('analysis/${nb}.Rmd', output_dir='$OUTDIR',
                                  output_file='${nb}.html', quiet=TRUE)" \
        > "$OUTDIR/${nb}.log" 2>&1
    st=$?
    check "${nb}.Rmd renders" "$st"
    [[ "$st" -ne 0 ]] && tail -15 "$OUTDIR/${nb}.log"
    [[ -s "$OUTDIR/${nb}.html" ]]; check "${nb}.html written" $?
done

step "5. Check the notebooks produced markers, not just the empty-marker guard"
for sample in CVB3_TT EVA71_TT; do
    n=$(Rscript -e "
      f <- file.path('$RES', '${sample}.markers_Pos.csv')
      m <- tryCatch(read.csv(f), error = function(e) data.frame())
      cat(nrow(m))" 2>/dev/null)
    [[ "${n:-0}" -gt 0 ]]
    check "${sample}: ${n:-0} marker rows (real path, not the guard)" $?
done

echo
echo "################################################################"
if [[ "$failures" -eq 0 ]]; then
    echo "# Seurat smoke test: all passed."
else
    echo "# Seurat smoke test: $failures check(s) FAILED."
fi
echo "################################################################"
[[ "$failures" -eq 0 ]]
