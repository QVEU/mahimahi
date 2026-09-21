#!/usr/bin/env python3
"""Concatenate per-sample strand-count tables into one analysis-ready table.

Replaces CollectFiles() from Scissors_Analysis_v4.ipynb, which dropped the Neg
and Pos columns (IF[,-c("Neg","Pos")]) so that a later cell plotting Neg/Pos
could not run. Nothing is dropped here.
"""

import argparse
import sys

import pandas as pd


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("tables", nargs="+", help="per-sample .tsv.gz files")
    p.add_argument("-o", "--output", required=True)
    args = p.parse_args(argv)

    frames, expected = [], None
    for path in args.tables:
        frame = pd.read_table(path)
        if frame.empty:
            print(f"  WARNING {path} is empty; skipping", file=sys.stderr)
            continue
        if expected is None:
            expected = list(frame.columns)
        elif list(frame.columns) != expected:
            sys.exit(f"merge_samples: column mismatch in {path}.\n"
                     f"  expected: {expected}\n  found:    {list(frame.columns)}")
        frames.append(frame)

    if not frames:
        sys.exit("merge_samples: every input table was empty.")

    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(args.output, sep="\t", index=False, compression="gzip")

    print(f"  merged {len(frames)} samples -> {len(merged):,} rows, "
          f"{merged['CBC'].nunique():,} cells, "
          f"{merged['ref_name'].nunique()} templates", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
