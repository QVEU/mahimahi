#!/usr/bin/env python3
"""Stage 4: doublet scoring, emitting barcode-keyed scores.

Usage:
    python3 scripts/04_scrublet.py <output_id> [--expected-doublet-rate 0.06]

Writes <cellranger_dir>/<output_id>/outs/filtered_feature_bc_matrix/
       <output_id>_Doublet_scores.tsv

Uses `scanpy.pp.scrublet`, not the standalone `scrublet` package. Scrublet was
absorbed into scanpy, and the standalone distribution has been frozen at 0.2.3
for years while scanpy stays maintained -- so depending on it directly means
depending on an unmaintained package for a step in the middle of the pipeline.

NOTE ON PARAMETERS. The original scrublet.Rmd referenced in the analysis
comments was never committed, so its settings could not be recovered. scanpy's
own default expected_doublet_rate is 0.05; this script keeps 0.06 and passes it
explicitly, so the value in force is visible rather than inherited from
whatever the library currently defaults to.

The per-sample doublet_scores cutoffs in analysis/config.R (0.45-0.58) were
tuned against the original standalone scrublet and are only meaningful on the
same scale, so re-pick them from the distributions this script produces rather
than carrying the old numbers over.

Output format (tab-separated, with header):
    barcode <TAB> doublet_score <TAB> predicted_doublet

The barcode column is the important part. The original files carried only two
unnamed columns, which forced the R side to attach scores to cells by row
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


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("output_id", help="sample directory name under the CellRanger dir")
    p.add_argument("--cellranger-dir", default=DEFAULT_CELLRANGER_DIR)
    p.add_argument("--matrix-dir",
                   help="read this matrix directory directly, instead of deriving "
                        "it from --cellranger-dir and the sample id")
    p.add_argument("--output", help="output TSV path (default: inside the matrix dir)")

    p.add_argument("--expected-doublet-rate", type=float, default=0.06,
                   help="10x multiplet rate for the number of cells loaded "
                        "(default: 0.06; scanpy's own default is 0.05)")
    p.add_argument("--sim-doublet-ratio", type=float, default=2.0)
    p.add_argument("--n-prin-comps", type=int, default=30)
    p.add_argument("--threshold", type=float, default=None,
                   help="override the automatic doublet-score threshold")
    p.add_argument("--random-state", type=int, default=0,
                   help="fixed so reruns reproduce the same scores")
    p.add_argument("--batch-key", default=None,
                   help="obs column to score within, if the matrix holds several runs")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    try:
        import numpy as np
        import scanpy as sc
    except ImportError as exc:
        sys.exit(f"Missing dependency: {exc}.\n"
                 "Install the environment with:\n"
                 "    conda env create -f workflow/envs/scissors.yaml && "
                 "conda activate scissors")

    # scanpy.pp.scrublet needs scikit-image to pick the doublet-score threshold
    # automatically, and it is NOT a hard dependency of scanpy -- a plain
    # `pip install scanpy` leaves it out. Check up front rather than failing
    # part-way through a long run on a large matrix.
    if args.threshold is None:
        try:
            import skimage  # noqa: F401
        except ImportError:
            sys.exit(
                "scanpy.pp.scrublet needs scikit-image to choose the doublet-score "
                "threshold automatically, and it is not installed.\n"
                "Either install it:\n"
                "    conda install -c conda-forge scikit-image\n"
                "(it is included in workflow/envs/scissors.yaml)\n"
                "or set the threshold yourself with --threshold.")

    matrix_dir = args.matrix_dir or os.path.join(
        args.cellranger_dir, args.output_id, "outs", "filtered_feature_bc_matrix")
    if not os.path.isdir(matrix_dir):
        sys.exit(f"Matrix directory not found: {matrix_dir}")

    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    print(f"{args.output_id}: {adata.n_obs} barcodes x {adata.n_vars} features",
          file=sys.stderr)

    if adata.n_obs < 30:
        sys.exit(f"{args.output_id}: only {adata.n_obs} barcodes. Scrublet simulates "
                 "doublets from the observed cells and needs a real population; "
                 "this matrix is too small to score.")

    sc.pp.scrublet(
        adata,
        batch_key=args.batch_key,
        sim_doublet_ratio=args.sim_doublet_ratio,
        expected_doublet_rate=args.expected_doublet_rate,
        n_prin_comps=args.n_prin_comps,
        threshold=args.threshold,
        random_state=args.random_state,
        verbose=True,
    )

    ## scanpy writes these two columns into obs.
    for column in ("doublet_score", "predicted_doublet"):
        if column not in adata.obs:
            sys.exit(f"scanpy.pp.scrublet did not produce obs['{column}']; got "
                     f"{list(adata.obs.columns)}. Check the scanpy version.")

    scores = adata.obs["doublet_score"].to_numpy()
    predicted = adata.obs["predicted_doublet"]

    out_path = args.output or os.path.join(
        matrix_dir, f"{args.output_id}_Doublet_scores.tsv")

    with open(out_path, "w") as handle:
        handle.write("barcode\tdoublet_score\tpredicted_doublet\n")
        for barcode, score, pred in zip(adata.obs_names, scores, predicted):
            # Write True/False explicitly rather than pandas' repr, so the R
            # side sees a stable spelling regardless of the column's dtype.
            flag = "NA" if pred is None else ("True" if bool(pred) else "False")
            handle.write(f"{barcode}\t{score:.6f}\t{flag}\n")

    print(f"Wrote {len(scores)} scores to {out_path}", file=sys.stderr)
    print(f"Score range: {float(np.nanmin(scores)):.4f} - "
          f"{float(np.nanmax(scores)):.4f}", file=sys.stderr)
    print(f"Predicted doublets: {int(predicted.sum())} / {len(scores)} "
          f"(expected rate {args.expected_doublet_rate})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
