#!/usr/bin/env python3
"""Report strand balance per template, and optionally assert the convention.

Mapped orientation is not the same thing as the strand of the original RNA:
the relationship depends on library chemistry (10x 3' vs 5') and on which mate
was aligned. Getting it backwards swaps Pos and Neg across every template and
produces output that looks entirely reasonable, because the numbers are all
still there -- just relabelled.

This turns that into something checkable. If the reference contains a sequence
whose sense orientation is known -- a host transcript, a spike-in, or a
reporter cassette -- name it with --sense-control. Host mRNA is overwhelmingly
positive-sense, so if the convention is right that template must come out
almost entirely 'Pos'. If it comes out mostly 'Neg', the convention is
inverted and every downstream ratio is upside down.

For a (+)-strand RNA virus, expect the viral template to be strongly Pos too,
with a Neg minority of roughly 0.5-5%: the negative strand is the replication
intermediate, not the dominant species. A near-50/50 split on a viral template
usually means the strand call is not working rather than that replication is
extraordinary.
"""

import argparse
import sys

import pandas as pd


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reads", help="per-read .tsv.gz from extract_reads.py")
    p.add_argument("-o", "--output", required=True, help="per-template TSV report")
    p.add_argument("--sample", required=True)
    p.add_argument("--sense-control",
                   help="reference name whose reads must be predominantly positive-sense")
    p.add_argument("--sense-control-min-frac", type=float, default=0.90,
                   help="minimum Pos fraction required of --sense-control (default: 0.90)")
    p.add_argument("--warn-balanced-frac", type=float, default=0.35,
                   help="warn when any template's Neg fraction exceeds this (default: 0.35)")
    args = p.parse_args(argv)

    reads = pd.read_table(args.reads, usecols=["CBC", "UMI", "ref_name", "strand"])
    if reads.empty:
        sys.exit("strand_qc: no reads to summarise.")

    distinct = reads.drop_duplicates(subset=["CBC", "ref_name", "strand", "UMI"])
    table = (distinct.groupby(["ref_name", "strand"], observed=True).size()
             .unstack(fill_value=0).reset_index())
    table.columns.name = None
    for strand in ("Pos", "Neg"):
        if strand not in table.columns:
            table[strand] = 0

    table["total_umis"] = table["Pos"] + table["Neg"]
    table["pos_frac"] = table["Pos"] / table["total_umis"]
    table["neg_frac"] = table["Neg"] / table["total_umis"]
    table["neg_over_pos"] = table["Neg"] / table["Pos"].where(table["Pos"] > 0)
    table["sample"] = args.sample
    table = table.sort_values("total_umis", ascending=False)

    table.to_csv(args.output, sep="\t", index=False)

    print(f"\n  strand balance, {args.sample} (distinct UMIs):", file=sys.stderr)
    print(f"  {'template':24s} {'Pos':>10s} {'Neg':>9s} {'neg_frac':>9s} {'Neg/Pos':>9s}",
          file=sys.stderr)
    for _, row in table.iterrows():
        ratio = "n/a" if pd.isna(row["neg_over_pos"]) else f"{row['neg_over_pos']:.4f}"
        print(f"  {str(row['ref_name']):24s} {int(row['Pos']):>10,} {int(row['Neg']):>9,} "
              f"{row['neg_frac']:>9.4f} {ratio:>9s}", file=sys.stderr)

    problems = []
    if args.sense_control:
        row = table[table["ref_name"] == args.sense_control]
        if row.empty:
            problems.append(
                f"--sense-control '{args.sense_control}' is not a reference in this "
                f"alignment. Present: {sorted(table['ref_name'].astype(str))}")
        else:
            frac = float(row["pos_frac"].iloc[0])
            if frac < args.sense_control_min_frac:
                problems.append(
                    f"sense control '{args.sense_control}' is only {frac:.1%} "
                    f"positive-sense (require >={args.sense_control_min_frac:.0%}).\n"
                    f"    The strand convention is probably inverted. Re-run with the "
                    f"other setting of strand_convention in config.yaml.")
            else:
                print(f"\n  sense control '{args.sense_control}': {frac:.1%} Pos -- "
                      f"convention looks right.", file=sys.stderr)

    for _, row in table.iterrows():
        if row["total_umis"] >= 100 and row["neg_frac"] > args.warn_balanced_frac:
            print(f"  WARNING {row['ref_name']} is {row['neg_frac']:.1%} negative-sense, "
                  f"which is high for a (+)-strand virus; check the strand call.",
                  file=sys.stderr)

    if problems:
        sys.exit("\nstrand_qc FAILED:\n  - " + "\n  - ".join(problems))
    return 0


if __name__ == "__main__":
    sys.exit(main())
