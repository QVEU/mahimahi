#!/usr/bin/env python3
"""Generate synthetic mahimahi inputs with known ground truth.

Builds a three-sequence template (two reporters plus a host-style sense
control) and reads whose per-cell, per-template, per-strand UMI counts are
known exactly, so the workflow's output can be checked against them.

The counts are chosen so that two templates in the SAME cell have clearly
different Rep_Index values. That is the case the original pivot_table(index=
"CBC") collapsed: it averaged across templates and reported one number for
both, which is why donor and acceptor replication could not be distinguished.
"""

import gzip
import io
import json
import shutil
import subprocess
import sys
import os
import random

random.seed(20240921)

OUT = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(OUT, "fixtures")
os.makedirs(FIXTURES, exist_ok=True)

BARCODE_LEN, UMI_LEN = 16, 12
READ_LEN = 120

TEMPLATES = ["eGFP", "mRuby3", "HostControl"]

# Two references with IDENTICAL sequence. Reads drawn from them map equally
# well to both, so minimap2 gives them MAPQ 0 and the min_mapq filter must
# drop every one. No Paralog row may appear in the output.
PARALOGS = ["ParalogA", "ParalogB"]

# (cell, template, strand) -> number of DISTINCT UMIs
#
# CELL_A carries both reporters at different replication levels:
#     eGFP    Neg/(Pos+Neg) =  4/44 = 0.0909
#     mRuby3  Neg/(Pos+Neg) =  5/15 = 0.3333
# A per-template-correct pipeline reports both. The original reported one
# blended number for each.
TRUTH = {
    ("CELL_AAAAAAAAAAAAA", "eGFP"):        {"Pos": 40, "Neg": 4},
    ("CELL_AAAAAAAAAAAAA", "mRuby3"):      {"Pos": 10, "Neg": 5},
    ("CELL_AAAAAAAAAAAAA", "HostControl"): {"Pos": 30, "Neg": 0},
    ("CELL_BBBBBBBBBBBBB", "eGFP"):        {"Pos": 20, "Neg": 1},
    ("CELL_BBBBBBBBBBBBB", "HostControl"): {"Pos": 25, "Neg": 0},
    ("CELL_CCCCCCCCCCCCC", "mRuby3"):      {"Pos": 8,  "Neg": 2},
}

# A population of cells with a known per-template (-)/(+) rate. Gives the
# fixtures enough cells, at a spread of depths, for a downstream slope fit of
# Neg ~ Pos to recover the rate the reads were built from.
POPULATION_RATES = {"eGFP": 0.05, "mRuby3": 0.20}
POPULATION_CELLS = 24

for _n in range(POPULATION_CELLS):
    _cell = f"CELL_POP{_n:03d}"
    for _template, _rate in POPULATION_RATES.items():
        _pos = 60 + 20 * _n              # 60..520, gives the fit real leverage
        _neg = int(round(_rate * _pos))
        TRUTH[(_cell, _template)] = {"Pos": _pos, "Neg": _neg}
    TRUTH[(_cell, "HostControl")] = {"Pos": 40, "Neg": 0}

# Reads per distinct UMI. >1 exercises deduplication: the pipeline must count
# distinct UMIs, not reads.
READS_PER_UMI = 3

# Reads that must NOT be counted. Each carries a real cell barcode, so a
# filter that leaks lets them through as extra counts or extra rows, which
# the verifier catches. The expected per-sample dispositions are written to
# expected.json and checked against each sample's extract_stats.json.
#
# FASTQ inputs (DRAGEN and raw 10x share one read set):
DECOY_UNMAPPED = 30    # random sequence: no alignment to any reference
DECOY_MULTIMAP = 20    # from the paralog pair: MAPQ 0
DECOY_NO_BARCODE = 10  # barcode block too short (DRAGEN: name field;
                       # raw 10x: R1, which embed_barcodes drops before alignment)
# Tagged BAM input: one count per excluded FLAG/MAPQ/tag case.
BAM_DECOYS_EACH = 10


def umi_for(index):
    """12 nt UMI from an integer, drawn from ONE pool shared by every cell,
    template and strand.

    The same UMI string therefore appears on both strands of a template, on
    several templates in a cell, and in several cells. Those are different
    molecules, and each must be counted. A pipeline that deduplicates on too
    narrow a key (UMI alone, or CBC+UMI, or omitting strand) undercounts here,
    and the exact comparison fails.
    """
    return "".join("ACGT"[(index >> (2 * i)) & 3] for i in range(UMI_LEN))[::-1]


def gzip_text(path):
    """Open a .gz for text writing with a zeroed mtime.

    gzip.open() stamps the current time into the gzip header, so regenerating
    identical reads yields a different file and a spurious diff every time.
    These fixtures are committed, so the bytes have to be reproducible. The
    reads themselves are already deterministic (fixed seed) -- only the header
    moved.
    """
    return io.TextIOWrapper(
        gzip.GzipFile(filename=path, mode="wb", compresslevel=6, mtime=0),
        encoding="utf-8")


def _cbc_for(cell):
    """Stable 16 nt barcode per cell name, padded with a non-informative base."""
    stem = cell.replace("CELL_", "").replace("POP", "P")
    return stem[:BARCODE_LEN].ljust(BARCODE_LEN, "A")


def random_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main():
    # --- template FASTA ------------------------------------------------
    refs = {name: random_seq(800) for name in TEMPLATES}
    paralog_seq = random_seq(800)
    for name in PARALOGS:
        refs[name] = paralog_seq
    template_path = os.path.join(FIXTURES, "template.fa")
    with open(template_path, "w") as fh:
        for name, seq in refs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # --- reads ----------------------------------------------------------
    # strand convention: reverse_is_positive. A read that maps REVERSE is
    # positive-sense, so a "Pos" read carries the reverse complement of the
    # reference fragment.
    records = []   # (read_id, cbc, umi, seq)
    counter = 0
    for (cell, template), strands in TRUTH.items():
        cbc = _cbc_for(cell)
        for strand, n_umis in strands.items():
            for u in range(n_umis):
                umi = umi_for(u)
                start = random.randint(0, len(refs[template]) - READ_LEN)
                fragment = refs[template][start:start + READ_LEN]
                seq = revcomp(fragment) if strand == "Pos" else fragment
                for _ in range(READS_PER_UMI):
                    counter += 1
                    records.append((f"SIM:1:FLOWCELL:1:1101:{counter}:{counter}",
                                    cbc, umi, seq))

    # Decoys. All use CELL_A's barcode and UMIs past any real one, so if one
    # leaks it shows up as an extra count or an extra row.
    decoy_cbc = _cbc_for("CELL_AAAAAAAAAAAAA")
    decoy_umi = iter(range(100_000, 200_000))
    for _ in range(DECOY_UNMAPPED):
        counter += 1
        records.append((f"SIM:1:FLOWCELL:1:1101:{counter}:{counter}",
                        decoy_cbc, umi_for(next(decoy_umi)), random_seq(READ_LEN)))
    for _ in range(DECOY_MULTIMAP):
        counter += 1
        start = random.randint(0, len(paralog_seq) - READ_LEN)
        records.append((f"SIM:1:FLOWCELL:1:1101:{counter}:{counter}",
                        decoy_cbc, umi_for(next(decoy_umi)),
                        paralog_seq[start:start + READ_LEN]))
    for _ in range(DECOY_NO_BARCODE):
        counter += 1
        start = random.randint(0, len(refs["eGFP"]) - READ_LEN)
        # Truncated barcode block: 16 nt CBC but only 4 nt of UMI.
        records.append((f"SIM:1:FLOWCELL:1:1101:{counter}:{counter}",
                        decoy_cbc, umi_for(next(decoy_umi))[:4],
                        refs["eGFP"][start:start + READ_LEN]))
    random.shuffle(records)

    qual = "I" * READ_LEN

    # --- (a) DRAGEN-style: CBC+UMI in colon field 7 of the read name ----
    with gzip_text(os.path.join(FIXTURES, "dragen_S1_R2_001.fastq.gz")) as fh:
        for rid, cbc, umi, seq in records:
            fh.write(f"@{rid}:{cbc}{umi}\n{seq}\n+\n{qual}\n")

    # --- (b) raw 10x: barcode in R1, cDNA in R2 -------------------------
    with gzip_text(os.path.join(FIXTURES, "raw10x_S2_R1_001.fastq.gz")) as r1, \
         gzip_text(os.path.join(FIXTURES, "raw10x_S2_R2_001.fastq.gz")) as r2:
        for rid, cbc, umi, seq in records:
            # barcode + UMI + polyT tail; a truncated UMI leaves R1 short.
            bc_seq = cbc + umi + ("T" * 10 if len(umi) == UMI_LEN else "")
            r1.write(f"@{rid} 1:N:0:1\n{bc_seq}\n+\n{'I' * len(bc_seq)}\n")
            r2.write(f"@{rid} 2:N:0:1\n{seq}\n+\n{qual}\n")

    # --- (c) Cell Ranger-style BAM with CB/UB tags ----------------------
    # Written as SAM here; the workflow's normalise_alignments sorts it.
    sam_path = os.path.join(FIXTURES, "tagged_S3.sam")
    bam_path = os.path.join(FIXTURES, "tagged_S3.bam")
    with open(sam_path, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        for name, seq in refs.items():
            fh.write(f"@SQ\tSN:{name}\tLN:{len(seq)}\n")
        for (cell, template), strands in TRUTH.items():
            cbc = _cbc_for(cell)
            for strand, n_umis in strands.items():
                for u in range(n_umis):
                    umi = umi_for(u)
                    for r in range(READS_PER_UMI):
                        counter += 1
                        flag = 16 if strand == "Pos" else 0
                        pos = random.randint(1, len(refs[template]) - READ_LEN)
                        frag = refs[template][pos - 1:pos - 1 + READ_LEN]
                        s = revcomp(frag) if strand == "Pos" else frag
                        fh.write(
                            f"SIM_T:{counter}\t{flag}\t{template}\t{pos}\t60\t"
                            f"{READ_LEN}M\t*\t0\t0\t{s}\t{qual}\t"
                            f"CB:Z:{cbc}-1\tUB:Z:{umi}\n")

        # Records that must be excluded: (flag, mapq, tags) per case. Each is
        # otherwise a perfectly good eGFP alignment in CELL_A.
        decoy_cbc = _cbc_for("CELL_AAAAAAAAAAAAA")
        bam_decoys = [
            (4,    0,  True),    # unmapped
            (256,  60, True),    # secondary
            (2048, 60, True),    # supplementary
            (0,    10, True),    # MAPQ below filters.min_mapq
            (0,    60, False),   # no UB tag
        ]
        decoy_umi = iter(range(100_000, 200_000))
        for flag, mapq, has_ub in bam_decoys:
            for _ in range(BAM_DECOYS_EACH):
                counter += 1
                pos = random.randint(1, len(refs["eGFP"]) - READ_LEN)
                s = refs["eGFP"][pos - 1:pos - 1 + READ_LEN]
                tags = f"CB:Z:{decoy_cbc}-1"
                if has_ub:
                    tags += f"\tUB:Z:{umi_for(next(decoy_umi))}"
                if flag & 4:
                    fh.write(f"SIM_T:{counter}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t"
                             f"{s}\t{qual}\t{tags}\n")
                else:
                    fh.write(f"SIM_T:{counter}\t{flag}\teGFP\t{pos}\t{mapq}\t"
                             f"{READ_LEN}M\t*\t0\t0\t{s}\t{qual}\t{tags}\n")

    # Convert to BAM and drop the SAM: 16 MB of text becomes a few hundred KB,
    # which is the difference between committable example data and not.
    if shutil.which("samtools"):
        # --no-PG: samtools otherwise records its own command line, including
        # the output path, in an @PG header line -- which makes the BAM differ
        # between runs in different directories and breaks --check.
        subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", bam_path, sam_path],
                       check=True)
        os.remove(sam_path)
    else:
        print("WARNING samtools not found; leaving tagged_S3.sam uncompressed. "
              "The committed fixture is a .bam, so run this where samtools "
              "is available before committing.", file=sys.stderr)

    # --- ground truth ---------------------------------------------------
    # Per-cell totals. UMI_count is distinct molecules, i.e. distinct
    # (template, strand, UMI); the shared UMI pool makes that differ from the
    # number of distinct UMI strings.
    cell_umis = {}
    for (cell, _template), strands in TRUTH.items():
        cbc = _cbc_for(cell)
        cell_umis[cbc] = cell_umis.get(cbc, 0) + sum(strands.values())

    kept = READS_PER_UMI * sum(cell_umis.values())
    dispositions = {
        "dragen_S1": {"kept": kept, "unmapped": DECOY_UNMAPPED,
                      "low_mapq": DECOY_MULTIMAP, "no_barcode": DECOY_NO_BARCODE},
        "raw10x_S2": {"kept": kept, "unmapped": DECOY_UNMAPPED,
                      "low_mapq": DECOY_MULTIMAP, "no_barcode": 0},
        "tagged_S3": {"kept": kept, "unmapped": BAM_DECOYS_EACH,
                      "secondary": BAM_DECOYS_EACH, "supplementary": BAM_DECOYS_EACH,
                      "low_mapq": BAM_DECOYS_EACH, "no_barcode": BAM_DECOYS_EACH},
    }

    expected = []
    for (cell, template), strands in TRUTH.items():
        cbc = _cbc_for(cell)
        pos, neg = strands.get("Pos", 0), strands.get("Neg", 0)
        expected.append({
            "CBC": cbc, "ref_name": template, "Pos": pos, "Neg": neg,
            "UMI_count": cell_umis[cbc],
            "CBC_readcount": READS_PER_UMI * cell_umis[cbc],
            "Rep_Index": (neg / (pos + neg)) if (pos + neg) else None,
            "Neg_PosRatio": neg / (pos + 1),
        })
    with open(os.path.join(FIXTURES, "expected.json"), "w") as fh:
        json.dump({"reads_per_umi": READS_PER_UMI,
                   "barcode_length": BARCODE_LEN, "umi_length": UMI_LEN,
                   "population_rates": POPULATION_RATES,
                   "population_cells": POPULATION_CELLS,
                   "dispositions": dispositions,
                   "expected": expected}, fh, indent=2)

    print(f"fixtures in {FIXTURES}")
    print(f"  templates: {TEMPLATES}")
    print(f"  {kept:,} countable reads, {kept // READS_PER_UMI:,} distinct molecules, "
          f"plus {len(records) - kept:,} FASTQ decoys")
    print(f"  {len(expected)} (cell, template) combinations")
    print(f"  population: {POPULATION_CELLS} cells at "
          + ", ".join(f"{k} (-)/(+)={v}" for k, v in POPULATION_RATES.items()))
    print("\n  ground truth (first 8 rows):")
    for e in expected[:8]:
        ri = "n/a" if e["Rep_Index"] is None else f"{e['Rep_Index']:.4f}"
        print(f"    {e['CBC'][:10]:10s} {e['ref_name']:12s} "
              f"Pos={e['Pos']:3d} Neg={e['Neg']:2d}  Rep_Index={ri}")


def check():
    """Regenerate into a temp dir and diff against the committed fixtures.

    The fixtures are committed, so they can silently fall out of step with the
    generator. This makes that a test failure rather than a surprise.
    """
    import filecmp
    import tempfile

    global FIXTURES
    committed = FIXTURES
    if not os.path.isdir(committed):
        sys.exit(f"--check: no committed fixtures at {committed}")
    expected = sorted(os.listdir(committed))

    with tempfile.TemporaryDirectory() as tmp:
        FIXTURES = tmp
        main()
        produced = sorted(os.listdir(tmp))

        if produced != expected:
            sys.exit("--check FAILED: file lists differ\n"
                     f"  committed: {expected}\n  regenerated: {produced}")

        differing = [f for f in expected
                     if not filecmp.cmp(os.path.join(committed, f),
                                        os.path.join(tmp, f), shallow=False)]

    if differing:
        sys.exit("--check FAILED: committed fixtures differ from regenerated:\n"
                 + "\n".join(f"  {f}" for f in differing)
                 + "\n\nRegenerate and commit:\n"
                   "  python3 tests/workflow/make_fixtures.py")
    print(f"--check OK: {len(expected)} committed fixtures are byte-identical "
          "to a fresh regeneration")


if __name__ == "__main__":
    if "--check" in sys.argv:
        check()
    else:
        main()
