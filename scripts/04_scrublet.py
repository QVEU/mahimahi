#!/usr/bin/env python3
"""Stage 4: doublet scoring with Scrublet, emitting barcode-keyed scores.

Usage:
    python scripts/04_scrublet.py <output_id> [--expected-doublet-rate 0.06]

Writes <cellranger_dir>/<output_id>/outs/filtered_feature_bc_matrix/
       <output_id>_Doublet_scores.tsv

NOTE: the original scrublet.Rmd referenced in the analysis comments was never
committed to the repository, so its parameters could not be recovered. This
script is a reimplementation using Scrublet's documented defaults. Before
reusing it on data you have already published, check that the scores it
produces line up with your existing *_Doublet_scores.tsv files -- the
per-sample doublet_scores cutoffs in the analysis scripts (0.45-0.58) were
tuned against the original output and are only meaningful on the same scale.

Output format (tab-separated, with header):
    barcode <TAB> doublet_score <TAB> predicted_doublet

The barcode column is the important change. The original files carried only
two unnamed columns, which forced the R side to attach scores to cells by row
position; any barcode dropped by CreateSeuratObject's min.features filter
shifted every subsequent score onto the wrong cell, silently. With barcodes
present, analysis/helpers.R joins by name instead.
"""

import argparse
import os
import sys

DEFAULT_CELLRANGER_DIR = os.environ.get(
    "SCISSORS_CELLRANGER_DIR",
    "/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_id", help="sample directory name under the CellRanger dir")
    parser.add_argument("--cellranger-dir", default=DEFAULT_CELLRANGER_DIR)
    parser.add_argument("--expected-doublet-rate", type=float, default=0.06,
                        help="10x multiplet rate for the loaded cell number (default 0.06)")
    parser.add_argument("--min-counts", type=int, default=2)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--n-prin-comps", type=int, default=30)
    parser.add_argument("--random-state", type=int, default=0,
                        help="fixed so reruns reproduce the same scores")
    args = parser.parse_args()

    try:
        import numpy as np
        import scanpy as sc
        import scrublet
    except ImportError as exc:
        print(f"Missing dependency: {exc}. Try: pip install scanpy scrublet", file=sys.stderr)
        return 1

    matrix_dir = os.path.join(
        args.cellranger_dir, args.output_id, "outs", "filtered_feature_bc_matrix"
    )
    if not os.path.isdir(matrix_dir):
        print(f"Matrix directory not found: {matrix_dir}", file=sys.stderr)
        return 1

    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
    print(f"{args.output_id}: {adata.n_obs} barcodes x {adata.n_vars} features")

    scrub = scrublet.Scrublet(
        adata.X,
        expected_doublet_rate=args.expected_doublet_rate,
        random_state=args.random_state,
    )
    scores, predicted = scrub.scrub_doublets(
        min_counts=args.min_counts,
        min_cells=args.min_cells,
        n_prin_comps=args.n_prin_comps,
    )

    if predicted is None:
        # Scrublet could not find a bimodal threshold; scores are still usable
        # with a manually chosen cutoff, which is how the analysis scripts work.
        print("Scrublet returned no automatic threshold; predicted_doublet set to NA.",
              file=sys.stderr)
        predicted_col = ["NA"] * len(scores)
    else:
        predicted_col = ["True" if p else "False" for p in predicted]

    out_path = os.path.join(matrix_dir, f"{args.output_id}_Doublet_scores.tsv")
    with open(out_path, "w") as handle:
        handle.write("barcode\tdoublet_score\tpredicted_doublet\n")
        for barcode, score, pred in zip(adata.obs_names, scores, predicted_col):
            handle.write(f"{barcode}\t{score:.6f}\t{pred}\n")

    print(f"Wrote {len(scores)} scores to {out_path}")
    print(f"Score range: {float(np.min(scores)):.4f} - {float(np.max(scores)):.4f}")
    if predicted is not None:
        print(f"Predicted doublets: {int(np.sum(predicted))} / {len(scores)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
