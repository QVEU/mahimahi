#!/usr/bin/env python3
"""Generate synthetic SCISSORS inputs with known ground truth.

Builds a three-sequence template (two reporters plus a host-style sense
control) and reads whose per-cell, per-template, per-strand UMI counts are
known exactly, so the workflow's output can be checked against them.

The counts are chosen so that two templates in the SAME cell have clearly
different Rep_Index values. That is the case the original pivot_table(index=
"CBC") collapsed: it averaged across templates and reported one number for
both, which is why donor and acceptor replication could not be distinguished.
"""

import gzip
import json
import os
import random

random.seed(20240921)

OUT = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(OUT, "fixtures")
os.makedirs(FIXTURES, exist_ok=True)

BARCODE_LEN, UMI_LEN = 16, 12
READ_LEN = 120

TEMPLATES = ["eGFP", "mRuby3", "HostControl"]

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

# A population of cells with a known per-template (-)/(+) rate, so the
# regression slope recovered by analysis/replication.R can be checked against
# the value the reads were built from. Kept deliberately noise-free: the point
# is to verify the plumbing from FASTQ through to a fitted slope, not to
# characterise the estimator (tests/test_replication_fit.R does that with
# noise).
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
                umi = f"{template[:2].upper()}{strand[0]}{u:04d}".ljust(UMI_LEN, "T")[:UMI_LEN]
                start = random.randint(0, len(refs[template]) - READ_LEN)
                fragment = refs[template][start:start + READ_LEN]
                seq = revcomp(fragment) if strand == "Pos" else fragment
                for _ in range(READS_PER_UMI):
                    counter += 1
                    records.append((f"SIM:1:FLOWCELL:1:1101:{counter}:{counter}",
                                    cbc, umi, seq))
    random.shuffle(records)

    qual = "I" * READ_LEN

    # --- (a) DRAGEN-style: CBC+UMI in colon field 7 of the read name ----
    with gzip.open(os.path.join(FIXTURES, "dragen_S1_R2_001.fastq.gz"), "wt") as fh:
        for rid, cbc, umi, seq in records:
            fh.write(f"@{rid}:{cbc}{umi}\n{seq}\n+\n{qual}\n")

    # --- (b) raw 10x: barcode in R1, cDNA in R2 -------------------------
    with gzip.open(os.path.join(FIXTURES, "raw10x_S2_R1_001.fastq.gz"), "wt") as r1, \
         gzip.open(os.path.join(FIXTURES, "raw10x_S2_R2_001.fastq.gz"), "wt") as r2:
        for rid, cbc, umi, seq in records:
            bc_seq = cbc + umi + "T" * 10          # barcode + UMI + polyT tail
            r1.write(f"@{rid} 1:N:0:1\n{bc_seq}\n+\n{'I' * len(bc_seq)}\n")
            r2.write(f"@{rid} 2:N:0:1\n{seq}\n+\n{qual}\n")

    # --- (c) Cell Ranger-style BAM with CB/UB tags ----------------------
    # Written as SAM here; the workflow's normalise_alignments sorts it.
    sam_path = os.path.join(FIXTURES, "tagged_S3.sam")
    with open(sam_path, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        for name, seq in refs.items():
            fh.write(f"@SQ\tSN:{name}\tLN:{len(seq)}\n")
        for (cell, template), strands in TRUTH.items():
            cbc = _cbc_for(cell)
            for strand, n_umis in strands.items():
                for u in range(n_umis):
                    umi = f"{template[:2].upper()}{strand[0]}{u:04d}".ljust(UMI_LEN, "T")[:UMI_LEN]
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

    # --- ground truth ---------------------------------------------------
    expected = []
    for (cell, template), strands in TRUTH.items():
        cbc = _cbc_for(cell)
        pos, neg = strands.get("Pos", 0), strands.get("Neg", 0)
        expected.append({
            "CBC": cbc, "ref_name": template, "Pos": pos, "Neg": neg,
            "Rep_Index": (neg / (pos + neg)) if (pos + neg) else None,
            "Neg_PosRatio": neg / (pos + 1),
        })
    with open(os.path.join(FIXTURES, "expected.json"), "w") as fh:
        json.dump({"reads_per_umi": READS_PER_UMI,
                   "barcode_length": BARCODE_LEN, "umi_length": UMI_LEN,
                   "population_rates": POPULATION_RATES,
                   "population_cells": POPULATION_CELLS,
                   "expected": expected}, fh, indent=2)

    print(f"fixtures in {FIXTURES}")
    print(f"  templates: {TEMPLATES}")
    print(f"  {len(records):,} reads, {len(records)//READS_PER_UMI:,} distinct UMIs")
    print(f"  {len(expected)} (cell, template) combinations")
    print(f"  population: {POPULATION_CELLS} cells at "
          + ", ".join(f"{k} (-)/(+)={v}" for k, v in POPULATION_RATES.items()))
    print("\n  ground truth (first 8 rows):")
    for e in expected[:8]:
        ri = "n/a" if e["Rep_Index"] is None else f"{e['Rep_Index']:.4f}"
        print(f"    {e['CBC'][:10]:10s} {e['ref_name']:12s} "
              f"Pos={e['Pos']:3d} Neg={e['Neg']:2d}  Rep_Index={ri}")


if __name__ == "__main__":
    main()
