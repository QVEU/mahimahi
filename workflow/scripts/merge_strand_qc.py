#!/usr/bin/env python3
"""Combine per-sample strand-balance reports into one summary table."""

import argparse
import sys

import pandas as pd


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reports", nargs="+", help="per-sample strand_qc.tsv files")
    p.add_argument("-o", "--output", required=True)
    args = p.parse_args(argv)

    frames = [pd.read_table(path) for path in args.reports]
    frames = [f for f in frames if not f.empty]
    if not frames:
        sys.exit("merge_strand_qc: every input report was empty.")

    summary = pd.concat(frames, ignore_index=True)
    summary.to_csv(args.output, sep="\t", index=False)

    print(f"  strand QC across {len(frames)} samples, "
          f"{len(summary)} template rows", file=sys.stderr)
    print(summary.to_string(index=False), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
