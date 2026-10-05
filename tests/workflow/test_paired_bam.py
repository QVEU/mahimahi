#!/usr/bin/env python3
"""Paired-end BAM: one molecule must count on one strand, not both.

A 10x 5' paired-end Cell Ranger BAM carries both mates of each read pair, with
CB/UB tags on each, and the mates map in opposite orientations. Strand is
defined by the cDNA read (R2), the same read the FASTQ route aligns. Counting
R1 as well would log every molecule once as Pos and once as Neg.

    python3 tests/workflow/test_paired_bam.py
"""

import gzip
import os
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
EXTRACT = os.path.join(HERE, "..", "..", "workflow", "scripts", "extract_reads.py")


def make_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "SQ": [{"SN": "PV", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        # Three molecules, each a mapped pair. Positive-sense RNA under the
        # reverse_is_positive convention: R2 reverse, R1 forward.
        for i, (umi, r2_reverse) in enumerate([("AAAAAAAAAAAA", True),
                                               ("CCCCCCCCCCCC", True),
                                               ("GGGGGGGGGGGG", False)]):
            for read1 in (True, False):
                a = pysam.AlignedSegment()
                a.query_name = f"pair{i}"
                a.query_sequence = "A" * 50
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                a.reference_id = 0
                a.reference_start = 100 + 10 * i
                a.cigarstring = "50M"
                a.mapping_quality = 60
                a.is_paired = True
                a.is_proper_pair = True
                a.is_read1 = read1
                a.is_read2 = not read1
                # Mates map in opposite orientations.
                a.is_reverse = r2_reverse if not read1 else not r2_reverse
                a.mate_is_reverse = not a.is_reverse
                a.next_reference_id = 0
                a.next_reference_start = a.reference_start
                a.set_tag("CB", "TTTTTTTTTTTTTTTT-1")
                a.set_tag("UB", umi)
                out.write(a)


def main():
    tmp = tempfile.mkdtemp(prefix="mahimahi_paired_")
    bam = os.path.join(tmp, "paired.bam")
    tsv = os.path.join(tmp, "reads.tsv.gz")
    make_bam(bam)
    proc = subprocess.run([sys.executable, EXTRACT, bam, "-o", tsv,
                           "--barcode-source", "tags"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stderr)
        sys.exit("FAILED: extract_reads.py exited non-zero")

    with gzip.open(tsv, "rt") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh][1:]
    strands = {}
    for read_id, cbc, umi, ref, strand, flag, mapq in rows:
        strands.setdefault(umi, set()).add(strand)

    failures = []
    if len(rows) != 3:
        failures.append(f"expected 3 rows (one per pair), got {len(rows)}")
    expect = {"AAAAAAAAAAAA": {"Pos"}, "CCCCCCCCCCCC": {"Pos"},
              "GGGGGGGGGGGG": {"Neg"}}
    for umi, want in expect.items():
        if strands.get(umi) != want:
            failures.append(f"UMI {umi}: strands {sorted(strands.get(umi, []))}, "
                            f"expected {sorted(want)}")

    if failures:
        print("FAILED:")
        for f in failures:
            print("  -", f)
        sys.exit(1)
    print("paired-end BAM: one strand per molecule, from R2 (3 pairs -> 2 Pos, 1 Neg)")


if __name__ == "__main__":
    main()
