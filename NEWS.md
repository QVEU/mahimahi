# Changes

## Narrowed to FASTQ → BAM → counts (October 2026)

mahimahi now does one thing: turn sequenced FASTQs (or an aligned BAM) into a
per-cell, per-template, per-strand UMI count table. The cell-analysis side has
been removed, because it was built from analysis scripts hardcoded to the
samples of earlier experiments and did not generalise. It is all in git
history (the last commit with it is `985a89c`) if anything needs to be
recovered.

Removed:

- `analysis/`: the Seurat analyses (`PV_mutants_integrated.R`,
  `CVB3_QC.Rmd`, `EVA71_QC.Rmd`), the replication-slope fitting and figures
  (`replication.R`, `mahimahi_replication.R`), and their per-sample config.
- `scripts/` and `config.sh`: the Cell Ranger `mkref` / `mkfastq` / `count` /
  `aggr` Slurm scripts, with their hardcoded Skyline paths and sample table,
  and the scrublet doublet step. Cell Ranger BAMs are still accepted as input
  (`barcode_source: tags`).
- `legacy/`: the original notebooks, and `docs/MIGRATION.md`, which mapped
  them to their replacements.
- The R unit tests, the Seurat smoke test and its stubs, and the CI job that
  ran them. The workflow's synthetic ground-truth tests are unchanged and
  remain the proof that counting works; the R slope check they used to run is
  gone with the R code.
- Seurat, R, scanpy and scikit-image from `workflow/envs/mahimahi.yaml`. The
  env is now seven packages and solves natively on Apple Silicon (the
  `r-mixtools` osx-arm64 gap no longer applies). CI runs the workflow tests on
  Linux and macOS.

Added:

- `outputs.keep_bam` in `config/config.yaml`: keep the minimap2 BAM and its
  index for FASTQ inputs instead of deleting it once counted.

## The strand-specific counting workflow replaces the original scripts

A Snakemake workflow replacing `mahimahi.sh`, `mahimahi_dragen.py`, and the
four divergent copies of the `SCISSORS()` tabulation
(`process_mahimahi.ipynb`, `process_mahimahi-Working020724.ipynb`,
`mahimahi.Rnb.ipynb`, `mahimahi_81423.ipynb`). Generalised to take FASTQ or an
already-aligned SAM/BAM/CRAM, mixed in one run, against any template FASTA.

**Per-template resolution was being destroyed by the pivot.**

Both Python `SCISSORS()` versions did
`pivot_table(index="CBC", columns="det_strand", values="UMI_strand_count")`.
pandas' default `aggfunc` is `"mean"`, and `UMI_strand_count` is repeated once
per read, so `Neg` and `Pos` became read-count-weighted averages across every
template in the cell. Worked example, one cell:

        truth    eGFP   Pos=1000 Neg=50   Rep_Index=0.0476
                 mRuby3 Pos= 100 Neg=20   Rep_Index=0.1667
        original both templates   Neg=41.43 Pos=918.18 Rep_Index=0.0432

`Rep_Index` and `Neg_PosRatio` came out identical for every template in a
cell, so they could not distinguish donor from acceptor replication. The
workflow keys the wide table on `(CBC, ref_name)`.

**`UMI_strand_count` grouping disagreed between versions.**
`process_mahimahi.ipynb` grouped by `("CBC","flag")`, omitting `ref_name`, so
even that column carried no per-template information.
`process_mahimahi-Working020724.ipynb` added `ref_name`.
`mahimahi.Rnb.ipynb` got it right via
`dcast(CBC+strand ~ ref_name, fun.aggregate = length(unique(UMI)))`.

**Barcode and UMI lengths are configuration.** `mahimahi_dragen.py` hardcoded
`[0:16]` and `[16:26]`. 10x 3' v2 is 16+10 but v3/v3.1 is 16+12, so on v3 the
UMI was truncated to 10 nt, shrinking the UMI space 16x. Under 1% undercount
at typical depth, but concentrated on the abundant `Pos` strand, so it inflates
`Neg/Pos` systematically.

**Alignment selection uses FLAG bits.** `flag.isin([0,16])` is equivalent to
"primary, unpaired, mapped" -- correct for single-end minimap2 output, but it
discards nearly everything in a paired or Cell Ranger BAM. `is_unmapped`,
`is_secondary`, `is_supplementary` and `is_reverse` generalise.

**Read names are validated before the run.** `split(":")[7]` raised
`IndexError` on any read name with fewer than 8 fields. The workflow samples
reads up front and fails with the offending name and the expected layout.

**One row per read.** `readDict[query_name] = ...` kept only the last record
for a read name, chosen by file order, whenever the aligner emitted more than
one.

**`mahimahi.py` (non-dragen) is not carried forward.** It calls `stop()`, which
does not exist in Python, and its `pd.DataFrame(readDict[i])` on a flat dict
raises `ValueError`. It could never have run.

**Strand convention is explicit and checkable.** Set `strand.sense_control` to
a sequence of known orientation and `rule strand_qc` fails the run if it does
not come out predominantly positive-sense.

**Output goes where it is told.** `to_csv(i.replace(".csv","_SCISSORS.csv"))`
had no directory component, so results landed in the interpreter's working
directory -- and `Mock_5h_S5_CBC.csv` appears in both the `EV71/` and `CVB/`
input directories, so the second run silently overwrote the first.

**Columns are selected by name**, not `.iloc[:, 12:]` / `.iloc[:, 3]`, which
depended on `pysam`'s `to_dict()` key order.

**R1 is not aligned.** `mahimahi.sh` ran minimap2 on the barcode read against
the viral template -- a 28 nt barcode has no meaningful viral alignment -- and
`mahimahi_dragen.py` then ignored that SAM entirely.

**`#!/bin/bash/`** in `mahimahi.sh` (trailing slash) is not a valid
interpreter path.

Tested end to end against synthetic reads with known per-cell, per-template,
per-strand UMI counts, over all three input routes:
`bash tests/workflow/run_workflow_tests.sh`.


## Renamed: SCISSORS to mahimahi

The pipeline's own names now say mahimahi. SCISSORS stays as the name of the
assay, and names that refer to the original scripts or to data on disk are
unchanged (`81123_SCISSORS_FreshAnalysis/`, `*_CBC_SCISSORS.csv`).

| Before | After |
|---|---|
| `workflow/envs/scissors.yaml`, conda env `scissors` | `workflow/envs/mahimahi.yaml`, conda env `mahimahi` |
| `results/scissors_counts.tsv.gz` | `results/mahimahi_counts.tsv.gz` |
| `SCISSORS_*` environment variables | `MAHIMAHI_*` (the old names are no longer read) |
