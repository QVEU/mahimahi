# Regenerating the strand counts on real data

The strand-counting workflow (`workflow/`) is the route for regenerating
`results/mahimahi_counts.tsv.gz`. This page records what has been verified on
the committed example data and lists what is still needed to run it on the
real samples.

## Verified on the example data

On a fresh clone, with only the workflow's tools installed:

```bash
conda create -n wf -c conda-forge -c bioconda \
    'python>=3.11' 'snakemake-minimal>=8.18' 'pysam>=0.22' 'pandas>=2.1' \
    'numpy>=1.24' 'minimap2>=2.26' 'samtools>=1.19' \
    'r-base>=4.4' 'r-ggplot2>=3.4' 'r-data.table>=1.14'
conda activate wf
bash tests/workflow/run_workflow_tests.sh
```

(The full `workflow/envs/mahimahi.yaml` works too; the list above is the
subset the workflow and its slope check need, and skips Seurat.)

Result, 2026-10-05, python 3.14.7, pysam 0.24.1, pandas 3.0.6, snakemake
9.27.0, minimap2 2.31, samtools 1.24:

- every per-cell, per-template `Pos`/`Neg`/`Rep_Index` value matches the
  synthetic truth for all three input routes (234 rows);
- eGFP and mRuby3 in the same cell come out distinct (0.0909 vs 0.3333);
- fitted slopes recover the generating rates (eGFP 0.04908, 95% CI
  0.04797–0.05018, true 0.05; mRuby3 0.19849, CI 0.19689–0.20010, true 0.20);
- an inverted strand convention is rejected by the sense control;
- a second run is a no-op, and a clean re-run gives a byte-identical
  decompressed `mahimahi_counts.tsv.gz`.

Without `Rscript` on PATH the slope check is skipped (or fails under
`--require-deps`); the count checks do not need R.

## What a real run needs

Nothing below is in the repository. `config/samples.tsv` and
`config/config.yaml` ship with placeholders (`/path/to/...`,
`resources/template.fa`), so the workflow cannot run until these are filled in.

### 1. One input per sample

For each sample, **one** of:

| You have | Sample sheet | `barcode_source` |
| --- | --- | --- |
| DRAGEN scRNA R2 FASTQ (barcode+UMI in the read name) | `input` = the R2 FASTQ | `read_name` (default) |
| Raw 10x FASTQ pair | `input` = R2, `mate_fastq` = R1 | `paired_fastq` |
| Cell Ranger `possorted_genome_bam.bam` | `input` = the BAM | `tags` |

The old per-read `*_CBC.csv` files under
`/hpcdata/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/81123_SCISSORS_FreshAnalysis/{CVB,EV71,PV_Rep}/`
are **not** usable inputs: they were produced by `mahimahi_dragen.py`, which
truncated v3 UMIs and kept only `flag in (0, 16)`. The workflow needs the
FASTQ or BAM upstream of them.

Samples the legacy notebooks analysed, as a starting list:

- **CVB3 / EVA71:** `Mock_5h_S5`, `CVB3_HH_S1`, `CVB3_TT_S2`, `EVA71_HH_S3`,
  `EVA71_TT_S4`
- **PV donor/acceptor:** `No-Acc_S1`, `No-Acc_S1i`, `C109S-Acc`, `D177A-Acc`,
  `dIRES-Acc`, `dIRESGAA-Acc`, `GAA-Acc`, `Y88P-Acc` (several appear under
  two `_S<n>` numbers, i.e. two runs; both runs are wanted, and `datalabel`
  will group them)
- **PV mutants:** `WT_GFP`, `WT_IRES_GFP`, `RFPMutPol`, `WT_GFP_RFPMutPol`,
  `RFP_488P`, `WT_GFP_RFP_488P`, `RFP_C109S`, `WT_GFP_RFP_C109S`,
  `RFP_D177A`, `WT_GFP_RFP_D177A`, `Del_IRES_mRuby3`,
  `Del_IRES_mRuby3_MutPol`, `WT_IRES_mRuby3_MutPol`,
  `WT_IRES_mRuby3_MutPol_WT_IRES_GFP`, `WT_IRES_GFP_Del_IRES_mRuby3`

Confirm which of these are wanted and where each one's FASTQ or BAM lives.

### 2. The template FASTA (FASTQ inputs only)

`template:` in `config/config.yaml`. One entry per template the counts should
distinguish: the viral genome or replicon for each experiment (CVB3, EVA71,
PV), plus `eGFP` / `mRuby3` (or whatever names the donor and acceptor carry)
for the co-infection runs. Each entry name becomes a `ref_name` in the output,
so name them as they should appear in figures.

The exact sequences used for the original runs are the ones to use. If the
CVB3, EVA71 and PV experiments used different templates, they need separate
config files (one `template` per run), e.g. `config/cvb3.yaml`,
`config/pv_donor_acceptor.yaml`, each passed with `--configfile`.

For Cell Ranger BAM inputs the FASTA is unused; the `ref_name` values come
from the BAM header (`PV`, `GFP`, `mRuby` in the `00_mkref.sh` reference).

### 3. A sense control

`strand.sense_control`: the name of one FASTA entry whose orientation is known
and which the reads cover, such as a host transcript (e.g. a stretch of
*ACTB* or *GAPDH* mRNA) or a reporter's coding sequence. Without it the
workflow still runs but cannot check the strand convention, and an inverted
convention swaps `Pos` and `Neg` everywhere while looking plausible. A
viral-only FASTA has no such entry, so add one.

### 4. Chemistry: 10x 5' v3 (answered 2026-10-05)

5' v3 is 16 nt barcode + 12 nt UMI, which is what `config/config.yaml`
already sets (`barcode_length: 16`, `umi_length: 12`). If any run used an
older 5' v1/v2 kit, that run needs `umi_length: 10`.

In a 5' library R2 is antisense to the RNA, so positive-sense viral RNA maps
reverse: `strand.convention: reverse_is_positive`, the default, is the right
setting. The sense control still checks it on the real data.

Cell Ranger BAMs from 5' paired-end runs carry both mates, which map in
opposite orientations; the workflow counts R2 only (R1 is reported as
`mate1` in `extract_stats.json`).

### 5. DRAGEN read-name layout (DRAGEN inputs only)

Which colon-separated field of the read name holds `<CBC><UMI>`; the default
is field 7. Copying the first read header of one R2 FASTQ is enough to check
this. The workflow validates the layout against the first 1000 reads and
stops with the offending read if it is wrong.

### 6. Barcode whitelist (recommended)

Per sample, the called-cell list: Cell Ranger's
`filtered_feature_bc_matrix/barcodes.tsv.gz` or DRAGEN's
`*.scRNA.barcodeSummary.tsv`. Goes in the `whitelist` column. This restricts
counting to real cells and is what makes the later join to the Seurat
metadata lossless.

## Running it

On Skyline, from the repository root:

```bash
conda env create -f workflow/envs/mahimahi.yaml && conda activate mahimahi
bash tests/workflow/run_workflow_tests.sh          # confirms the install
snakemake --cores 8 --configfile config/<run>.yaml -n   # check the plan
snakemake --cores 8 --configfile config/<run>.yaml
```

Then check `results/<sample>/strand_qc.tsv` (sense control `pos_frac` >= 0.9;
viral `neg_over_pos` roughly 0.005–0.05) and
`results/<sample>/extract_stats.json` (`kept` close to `total`) before using
the counts.
