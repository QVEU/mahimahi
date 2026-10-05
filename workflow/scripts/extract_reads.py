#!/usr/bin/env python3
"""Extract per-read cell barcode, UMI, template and strand from a SAM/BAM.

Replaces mahimahi_dragen.py. Differences that matter:

- Streams. The original built a dict of every passing read in memory and then
  converted it to a DataFrame; this writes rows as it reads them.
- Keeps one row per read. The original used ``readDict[query_name] = ...``,
  so when an aligner emitted more than one record for a read name only the
  last one survived, chosen by file order.
- Selects alignments by FLAG bits, not by ``flag in (0, 16)``. The original
  test is equivalent to "primary, unpaired, mapped" and happens to be right
  for single-end minimap2 output, but discards essentially everything in a
  paired or Cell Ranger BAM. Bit tests generalise.
- Barcode and UMI lengths are arguments, not slice literals. The original
  hardcoded ``[0:16]`` and ``[16:26]``, which silently truncates a 12 nt
  10x v3 UMI to 10 nt.
- Validates the read-name layout on a sample of reads and fails with a clear
  message instead of raising IndexError partway through a file.
- Does not emit seq/qual. The original wrote both to CSV, which dominated
  output size.
"""

import argparse
import gzip
import sys
from collections import Counter

import pysam


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("alignments", help="input .sam/.bam/.cram")
    p.add_argument("-o", "--output", required=True, help="output .tsv.gz")
    p.add_argument("--stats", help="write a JSON summary of read dispositions here")

    p.add_argument("--barcode-source", choices=["read_name", "tags"], default="read_name",
                   help="where the cell barcode and UMI live (default: read_name)")
    p.add_argument("--name-format", choices=["dragen", "suffix"], default="dragen",
                   help="read_name layout: 'dragen' = colon-delimited field holding "
                        "CBC+UMI concatenated; 'suffix' = name ends _<CBC>_<UMI> "
                        "(what embed_barcodes.py writes)")
    p.add_argument("--name-field", type=int, default=7,
                   help="0-based colon-delimited field for --name-format dragen (default: 7)")
    p.add_argument("--barcode-length", type=int, default=16)
    p.add_argument("--umi-length", type=int, default=12,
                   help="10x 3' v2 is 10, v3/v3.1 is 12 (default: 12)")
    p.add_argument("--cb-tag", default="CB", help="corrected-barcode tag (default: CB)")
    p.add_argument("--ub-tag", default="UB", help="corrected-UMI tag (default: UB)")

    p.add_argument("--min-mapq", type=int, default=55,
                   help="minimum MAPQ (default: 55; minimap2 caps at 60)")
    p.add_argument("--keep-secondary", action="store_true")
    p.add_argument("--keep-supplementary", action="store_true")
    p.add_argument("--reverse-is-positive", dest="reverse_is_positive",
                   action="store_true", default=True,
                   help="reverse-mapped reads are positive-sense viral RNA "
                        "(the convention every original script uses)")
    p.add_argument("--forward-is-positive", dest="reverse_is_positive",
                   action="store_false",
                   help="invert the strand convention")
    p.add_argument("--validate-reads", type=int, default=1000,
                   help="reads to inspect up front for name-layout sanity (default: 1000)")
    return p.parse_args(argv)


def strip_cell_suffix(barcode):
    """Drop a 10x '-1' style suffix so barcodes join to Cell Ranger's."""
    return barcode.split("-", 1)[0]


def make_barcode_reader(args):
    """Return fn(AlignedSegment) -> (cbc, umi) or (None, None)."""
    if args.barcode_source == "tags":
        cb, ub = args.cb_tag, args.ub_tag

        def from_tags(rec):
            try:
                return strip_cell_suffix(rec.get_tag(cb)), rec.get_tag(ub)
            except KeyError:
                return None, None
        return from_tags

    if args.name_format == "suffix":
        def from_suffix(rec):
            parts = rec.query_name.rsplit("_", 2)
            if len(parts) != 3:
                return None, None
            return parts[1], parts[2]
        return from_suffix

    field, bc_len, umi_len = args.name_field, args.barcode_length, args.umi_length
    need = bc_len + umi_len

    def from_dragen(rec):
        parts = rec.query_name.split(":")
        if len(parts) <= field:
            return None, None
        block = parts[field]
        if len(block) < need:
            return None, None
        return block[:bc_len], block[bc_len:need]
    return from_dragen


def validate_name_layout(path, args, get_barcode):
    """Inspect the first N reads so a bad layout fails immediately and loudly."""
    if args.barcode_source == "tags" or args.validate_reads <= 0:
        return
    seen = ok = 0
    example = None
    with pysam.AlignmentFile(path, "r" if path.endswith(".sam") else "rb",
                             check_sq=False) as fh:
        for rec in fh:
            if example is None:
                example = rec.query_name
            seen += 1
            cbc, _ = get_barcode(rec)
            if cbc is not None:
                ok += 1
            if seen >= args.validate_reads:
                break
    if seen == 0:
        return
    if ok == 0:
        need = args.barcode_length + args.umi_length
        sys.exit(
            f"extract_reads: could not parse a barcode from any of the first {seen} reads.\n"
            f"  example read name: {example!r}\n"
            f"  --name-format {args.name_format}"
            + (f" --name-field {args.name_field}" if args.name_format == "dragen" else "")
            + f" expects at least {need} characters "
              f"({args.barcode_length} barcode + {args.umi_length} UMI).\n"
            f"  If the barcode is in a BAM tag instead, use --barcode-source tags."
        )
    if ok < seen * 0.5:
        print(f"extract_reads: WARNING only {ok}/{seen} sampled reads yielded a barcode; "
              f"check --name-format/--name-field.", file=sys.stderr)


def main(argv=None):
    args = parse_args(argv)
    get_barcode = make_barcode_reader(args)
    validate_name_layout(args.alignments, args, get_barcode)

    positive = "Pos"
    negative = "Neg"
    counts = Counter()

    mode = "r" if args.alignments.endswith(".sam") else "rb"
    with pysam.AlignmentFile(args.alignments, mode, check_sq=False) as fh, \
         gzip.open(args.output, "wt", compresslevel=6) as out:
        out.write("read_id\tCBC\tUMI\tref_name\tstrand\tflag\tmapq\n")
        for rec in fh:
            counts["total"] += 1
            if rec.is_unmapped:
                counts["unmapped"] += 1
                continue
            if rec.is_secondary and not args.keep_secondary:
                counts["secondary"] += 1
                continue
            if rec.is_supplementary and not args.keep_supplementary:
                counts["supplementary"] += 1
                continue
            if rec.mapping_quality < args.min_mapq:
                counts["low_mapq"] += 1
                continue

            cbc, umi = get_barcode(rec)
            if cbc is None or umi is None:
                counts["no_barcode"] += 1
                continue

            # is_reverse reads FLAG bit 0x10 regardless of the other bits.
            strand = (positive if rec.is_reverse else negative) \
                if args.reverse_is_positive else \
                (negative if rec.is_reverse else positive)

            counts["kept"] += 1
            counts[f"kept_{strand}"] += 1
            out.write(f"{rec.query_name}\t{cbc}\t{umi}\t{rec.reference_name}\t"
                      f"{strand}\t{rec.flag}\t{rec.mapping_quality}\n")

    for key in ("total", "kept", "kept_Pos", "kept_Neg", "unmapped", "secondary",
                "supplementary", "low_mapq", "no_barcode"):
        print(f"  {key:15s} {counts[key]:>12,}", file=sys.stderr)

    if counts["kept"] == 0:
        sys.exit("extract_reads: no reads passed filtering; refusing to write an empty table.")

    if args.stats:
        import json
        with open(args.stats, "w") as fh:
            json.dump({
                "alignments": args.alignments,
                "strand_convention": ("reverse_is_positive" if args.reverse_is_positive
                                      else "forward_is_positive"),
                "barcode_source": args.barcode_source,
                "barcode_length": args.barcode_length,
                "umi_length": args.umi_length,
                "min_mapq": args.min_mapq,
                "counts": dict(counts),
            }, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
