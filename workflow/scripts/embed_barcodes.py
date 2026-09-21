#!/usr/bin/env python3
"""Embed the 10x cell barcode and UMI from R1 into the R2 read name.

For a raw 10x library the barcode is the first bases of R1, not part of the
read name, so it has to travel with R2 into the aligner. This appends
_<CBC>_<UMI> to each R2 name; extract_reads.py --name-format suffix reads it
back.

Both FASTQs are streamed in lockstep, which is safe because Illumina writes
mates in the same order, and keeps memory constant regardless of library size.
Order is verified per record rather than assumed: if the names diverge the run
stops instead of silently pairing the wrong barcode to the wrong read.

(The original mahimahi.sh aligned R1 against the viral template with minimap2
to get at its sequence. A 28 nt barcode read has no meaningful alignment to a
viral genome, and the script that consumed the output ignored the R1 SAM
entirely, so that alignment was wasted wall-clock.)
"""

import argparse
import gzip
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def writer(path):
    return gzip.open(path, "wt", compresslevel=6) if path.endswith(".gz") else open(path, "w")


def fastq_records(handle, path):
    """Yield (name_line, seq, plus, qual). Validates the 4-line structure."""
    while True:
        name = handle.readline()
        if not name:
            return
        seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
        if not qual:
            sys.exit(f"embed_barcodes: {path} ends mid-record; file is truncated.")
        if not name.startswith("@"):
            sys.exit(f"embed_barcodes: {path} is not valid FASTQ "
                     f"(expected '@' at record start, got {name[:20]!r}).")
        yield name.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def read_id(name_line):
    """'@A00:1:... 2:N:0:IDX' -> 'A00:1:...' (drop '@' and the mate suffix)."""
    return name_line[1:].split()[0]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--r1", required=True, help="barcode read (CBC+UMI)")
    p.add_argument("--r2", required=True, help="cDNA read")
    p.add_argument("-o", "--output", required=True, help="output R2 .fastq.gz")
    p.add_argument("--barcode-length", type=int, default=16)
    p.add_argument("--umi-length", type=int, default=12)
    args = p.parse_args(argv)

    need = args.barcode_length + args.umi_length
    written = skipped_short = 0

    with opener(args.r1) as h1, opener(args.r2) as h2, writer(args.output) as out:
        for (n1, s1, _, _), (n2, s2, p2, q2) in zip(
                fastq_records(h1, args.r1), fastq_records(h2, args.r2)):
            id1, id2 = read_id(n1), read_id(n2)
            if id1 != id2:
                sys.exit(f"embed_barcodes: R1/R2 out of order at record {written + skipped_short + 1}.\n"
                         f"  R1: {id1}\n  R2: {id2}\n"
                         f"Re-sort or re-demultiplex; pairing by position is not safe here.")
            if len(s1) < need:
                skipped_short += 1
                continue
            cbc = s1[:args.barcode_length]
            umi = s1[args.barcode_length:need]
            out.write(f"@{id2}_{cbc}_{umi}\n{s2}\n{p2}\n{q2}\n")
            written += 1

    print(f"  embedded {written:,} barcodes"
          + (f"; skipped {skipped_short:,} R1 reads shorter than {need} nt"
             if skipped_short else ""), file=sys.stderr)
    if written == 0:
        sys.exit("embed_barcodes: no records written. Check --barcode-length/--umi-length "
                 "against the actual R1 read length.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
