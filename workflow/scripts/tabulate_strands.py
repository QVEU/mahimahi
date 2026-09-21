#!/usr/bin/env python3
"""Collapse per-read records into strand-specific UMI counts per cell and template.

This replaces the SCISSORS() function from process_mahimahi*.ipynb. Four
behaviours are deliberately different, because the original produced numbers
that could not be interpreted per template:

1. UMIs are deduplicated within (CBC, ref_name, strand).

   process_mahimahi.ipynb grouped by ("CBC", "flag") only, so every template
   in a cell shared one count. process_mahimahi-Working020724.ipynb added
   ref_name, which is correct.

2. The wide table is keyed on (CBC, ref_name).

   Both originals did pivot_table(index="CBC", ...). pandas' default aggfunc
   is "mean", so Neg and Pos became averages across every template in the
   cell -- and because UMI_strand_count was repeated once per read, a
   read-count-weighted average at that. Worked example, one cell:

       truth    eGFP   Pos=1000 Neg=50   Rep_Index=0.0476
                mRuby3 Pos= 100 Neg=20   Rep_Index=0.1667
       original both templates   Neg=41.43 Pos=918.18 Rep_Index=0.0432

   Rep_Index came out identical for both templates, so it could not
   distinguish donor from acceptor replication.

3. Output columns are selected by name.

   The original used .iloc[:, 12:] and .iloc[:, 3], which land on the right
   columns only for one particular pysam to_dict() key order.

4. Output goes where it is told.

   The original called to_csv(i.replace(...)) with no directory, writing into
   the interpreter's working directory. Two input directories each containing
   a Mock_5h_S5_CBC.csv silently overwrote one another.
"""

import argparse
import sys

import pandas as pd

LONG_COLUMNS = ["CBC", "ref_name", "strand", "sample", "datalabel",
                "UMI_strand_count", "CBC_readcount", "UMI_count"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reads", help="per-read .tsv.gz from extract_reads.py")
    p.add_argument("--sample", required=True, help="sample identifier")
    p.add_argument("--datalabel", help="genotype label shared across sequencing runs "
                                       "(default: sample with a trailing _S<n> removed)")
    p.add_argument("--long-output", required=True, help="long-format .tsv.gz")
    p.add_argument("--wide-output", required=True, help="wide-format .tsv.gz")
    p.add_argument("--min-cell-umis", type=int, default=0,
                   help="drop cells with fewer than this many total UMIs (default: 0). "
                        "The originals filtered on CBC_readcount > 100 here and then "
                        "again on UMI_count > 100 in R; keep filtering in one place.")
    p.add_argument("--whitelist", help="optional file of accepted cell barcodes, one per "
                                       "line, or a DRAGEN barcodeSummary.tsv")
    return p.parse_args(argv)


def default_datalabel(sample):
    """'WT_GFP_S1' -> 'WT_GFP'. Mirrors the originals' i.split('_S')[0], but
    only strips a genuine _S<digits> tail so a sample like 'WT_GFP_Series2'
    is left alone."""
    import re
    return re.sub(r"_S\d+$", "", sample)


def load_whitelist(path):
    """Accept a bare barcode list or a DRAGEN barcodeSummary.tsv."""
    head = pd.read_table(path, nrows=5)
    if "Barcode" in head.columns:
        table = pd.read_table(path)
        if "Filter" in table.columns:
            table = table[table["Filter"] == "PASS"]
        barcodes = table["Barcode"]
    else:
        barcodes = pd.read_table(path, header=None).iloc[:, 0]
    return set(barcodes.astype(str).str.split("-").str[0])


def main(argv=None):
    args = parse_args(argv)
    datalabel = args.datalabel or default_datalabel(args.sample)

    reads = pd.read_table(args.reads, usecols=["CBC", "UMI", "ref_name", "strand"],
                          dtype={"CBC": "string", "UMI": "string",
                                 "ref_name": "category", "strand": "category"})
    if reads.empty:
        sys.exit(f"tabulate_strands: {args.reads} has no rows.")

    if args.whitelist:
        allowed = load_whitelist(args.whitelist)
        before = len(reads)
        reads = reads[reads["CBC"].isin(allowed)]
        print(f"  whitelist: {before:,} -> {len(reads):,} reads "
              f"({len(allowed):,} accepted barcodes)", file=sys.stderr)
        if reads.empty:
            sys.exit("tabulate_strands: no reads left after whitelist filtering. "
                     "Check that the whitelist barcodes match the read barcodes.")

    # Reads per cell, before any UMI collapsing.
    cbc_readcount = reads.groupby("CBC", observed=True).size().rename("CBC_readcount")

    # Deduplicate to distinct UMIs. drop_duplicates then size is equivalent to
    # nunique but avoids materialising per-group hash sets.
    distinct = reads.drop_duplicates(subset=["CBC", "ref_name", "strand", "UMI"])

    # Total distinct UMIs per cell, across templates and strands. This is the
    # cell's viral UMI content.
    umi_count = distinct.groupby("CBC", observed=True)["UMI"].nunique().rename("UMI_count")

    # The count that matters: distinct UMIs per cell, per template, per strand.
    long = (distinct
            .groupby(["CBC", "ref_name", "strand"], observed=True)
            .size()
            .rename("UMI_strand_count")
            .reset_index())

    long = long.merge(cbc_readcount, on="CBC").merge(umi_count, on="CBC")
    long["sample"] = args.sample
    long["datalabel"] = datalabel

    # Wide: one row per (CBC, ref_name). Both strands become columns. Missing
    # strands are true zeros -- a cell with no negative-strand UMIs for a
    # template has Neg == 0, not NaN.
    wide = (long
            .pivot_table(index=["CBC", "ref_name"], columns="strand",
                         values="UMI_strand_count", aggfunc="sum",
                         fill_value=0, observed=True)
            .reset_index())
    wide.columns.name = None
    for strand in ("Pos", "Neg"):
        if strand not in wide.columns:
            wide[strand] = 0

    wide = wide.merge(long[["CBC", "CBC_readcount", "UMI_count"]].drop_duplicates(),
                      on="CBC", how="left")

    # Neg_PosRatio keeps the originals' +1 in the denominator so a cell with no
    # positive-strand UMIs does not divide by zero. Rep_Index is the fraction
    # of that template's UMIs that are negative-sense, which is the quantity
    # the plot axes call "(-)strand/total(vRNA)".
    wide["Neg_PosRatio"] = wide["Neg"] / (wide["Pos"] + 1)
    total = wide["Pos"] + wide["Neg"]
    wide["Rep_Index"] = (wide["Neg"] / total).where(total > 0, other=pd.NA)

    wide["sample"] = args.sample
    wide["datalabel"] = datalabel

    if args.min_cell_umis > 0:
        keep = wide["UMI_count"] >= args.min_cell_umis
        cells_before = wide["CBC"].nunique()
        wide = wide[keep]
        long = long[long["CBC"].isin(set(wide["CBC"]))]
        print(f"  min-cell-umis {args.min_cell_umis}: "
              f"{cells_before:,} -> {wide['CBC'].nunique():,} cells", file=sys.stderr)

    long = long[LONG_COLUMNS]
    wide = wide[["CBC", "ref_name", "sample", "datalabel", "CBC_readcount",
                 "UMI_count", "Neg", "Pos", "Neg_PosRatio", "Rep_Index"]]

    long.to_csv(args.long_output, sep="\t", index=False, compression="gzip")
    wide.to_csv(args.wide_output, sep="\t", index=False, compression="gzip")

    print(f"  {args.sample}: {wide['CBC'].nunique():,} cells x "
          f"{wide['ref_name'].nunique()} templates -> {len(wide):,} rows",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
