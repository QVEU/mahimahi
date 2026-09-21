# SCISSORS

**S**ingle **C**ell **I**solated **S**trand **S**pecific **O**bservation of
**R**eplication **S**tate

Developed by the Quantitative Virology and Evolution Unit (QVEU) at NIAID. The
pipeline uses 10x Genomics single-cell RNA-seq to simultaneously quantify
enterovirus replication and the host transcriptional response within the same
cells, by counting viral and reporter reads as ordinary features against a
custom reference.

This repository is a restructured version of `QVEU/SCISSORS`. The analysis is
unchanged in intent; the reorganization fixes a set of bugs and removes
hardcoded paths and indices. See [NEWS.md](NEWS.md) for the full list of
changes and [Known open questions](#known-open-questions) for the two things
that need a decision from someone who ran the original experiments.

Shell scripts target the Skyline HPC cluster and its Slurm scheduler.

---

## Pipeline overview

```
                    scripts/00_mkref.sh
              custom reference: GRCh38-2020-A
                  + PV genome + GFP + mRuby3
                              |
                    scripts/01_mkfastq.sh
                  BCL run dir -> per-sample FASTQs
                              |
                    scripts/02_count.sh
                  FASTQs -> filtered_feature_bc_matrix
                    (one job per sample, array-capable)
                              |
              +---------------+---------------+
              |                               |
   scripts/03_aggregate.sh          scripts/04_scrublet.py
   (optional; Seurat path            doublet scores per cell,
    does not read this)              keyed by barcode
                                              |
              +-------------------------------+
              |                               |
  analysis/CVB3_QC.Rmd            analysis/PV_mutants_integrated.R
  analysis/EVA71_QC.Rmd           14 PV samples: infected-cell calling,
  per-sample QC and clustering     QC, merge, cell-cycle regression,
                                   clustering, markers
```

## Repository layout

| Path | Purpose |
| --- | --- |
| `config.sh` | Cluster paths, module versions, Slurm-derived core/memory helpers. Sourced by every shell script. |
| `scripts/00_mkref.sh` | Build the custom transcriptome. Run once. |
| `scripts/01_mkfastq.sh` | Demultiplex a BCL run directory. |
| `scripts/02_count.sh` | Count one sample, or all of them as a Slurm array job. |
| `scripts/03_aggregate.sh` | Optional `cellranger aggr` across samples. |
| `scripts/04_scrublet.py` | Doublet scoring, emitting barcode-keyed TSVs. |
| `scripts/samples.tsv` | `output_id` / `fastq_sample` / `fastq_dir` per sample. |
| `analysis/config.R` | Share-root resolution, sample lists, per-sample QC thresholds. |
| `analysis/helpers.R` | Doublet joining, infected-cell calling, QC filtering. |
| `analysis/CVB3_QC.Rmd` | Coxsackievirus B3 per-sample QC and clustering. |
| `analysis/EVA71_QC.Rmd` | Enterovirus A71 per-sample QC and clustering. |
| `analysis/PV_mutants_integrated.R` | Integrated poliovirus mutant analysis. |
| `tests/` | Offline regression tests for the correctness fixes. |

## Running it

### 0. Build the reference (once)

```bash
sbatch scripts/00_mkref.sh
```

Combines human GRCh38-2020-A with the poliovirus genome and the GFP and mRuby3
reporter sequences, so `PV`, `GFP` and `mRuby` become countable features.
Refuses to overwrite an existing reference, because every downstream count
matrix is tied to it.

### 1. Demultiplex

```bash
sbatch scripts/01_mkfastq.sh /path/to/bcl/run/dir QVEU0056
```

Needs `SampleSheet.csv` at the path in `config.sh`. That file is an Illumina
sample sheet as `cellranger mkfastq` expects it — a `[Data]` section with
`Lane`, `Sample`, and `Index` columns, where `Index` is a 10x index set name
such as `SI-GA-A1`. It is **not** in this repository; it lives next to the
data.

Record the resulting `outs/fastq_path` directory in `scripts/samples.tsv`.

### 2. Count

One sample:

```bash
sbatch scripts/02_count.sh Mock_5h_PV
```

All samples in `scripts/samples.tsv` as an array job:

```bash
sbatch --array=1-$(grep -vc '^#' scripts/samples.tsv) scripts/02_count.sh
```

The script refuses to start if the sample's `fastq_dir` is still `TODO`, if the
transcriptome is missing, or if no FASTQs match the sample prefix — Cell Ranger
will otherwise run happily to completion on zero reads and produce an empty
matrix that is only noticed much later in Seurat.

### 3. Doublet scores

```bash
python3 scripts/04_scrublet.py Mock_5h_PV
```

Writes `<sample>_Doublet_scores.tsv` into that sample's
`filtered_feature_bc_matrix/` with a `barcode` column. See
[Doublet scores](#doublet-scores-read-this-before-reusing-old-tsvs) below.

### 4. Aggregation CSV (optional)

`scripts/03_aggregate.sh` takes a CSV you write yourself:

```csv
sample_id,molecule_h5
Mock_5h_PV,/path/to/CellRanger/Mock_5h_PV/outs/molecule_info.h5
WT_GFP_PV,/path/to/CellRanger/WT_GFP_PV/outs/molecule_info.h5
```

```bash
sbatch scripts/03_aggregate.sh /path/to/aggr.csv mutants
```

The Seurat analyses do **not** read `aggr` output — they read the per-sample
matrices and merge in R. Run this only if you want Cell Ranger's own
depth-normalized matrix.

### 5. Analysis

```bash
Rscript analysis/PV_mutants_integrated.R
```

and knit `analysis/CVB3_QC.Rmd` / `analysis/EVA71_QC.Rmd` in RStudio.

Paths resolve automatically: `analysis/config.R` looks for the lab share at
`/data/lvd_qve` (Skyline), `~/lab_share`, or `/Volumes/lvd_qve` (macOS), in
that order. Override with `SCISSORS_SHARE_ROOT=/some/path`. Shell paths are
overridable with `SCISSORS_PROJECT_ROOT` and `SCISSORS_REF_BUILD_DIR`.

### Environment

The R analyses target **Seurat v4** and check this at startup. They will not
run unchanged on v5, where layers change `merge()` and RNA assay behaviour.

```r
remotes::install_version("SeuratObject", "4.1.4",
  repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0",
  repos = c("https://satijalab.r-universe.dev", getOption("repos")))
install.packages(c("clustree", "tidyverse", "ggrepel", "mixtools", "rprojroot"))
```

Cluster modules: `cellranger/7.2.0-dntehee`, `bcl2fastq2/2.20.0.422-orocbiu`.

## Doublet scores: read this before reusing old TSVs

The doublet score files produced by the original workflow had **two unnamed
columns** — score and prediction, no barcodes. The R code therefore attached
them to cells by row position:

```r
doublets <- read.table(path, header = FALSE)
obj <- AddMetaData(obj, metadata = doublets$doublet_scores, ...)
```

`AddMetaData` given an unnamed vector assigns it in object order. But the
Seurat objects were built with `min.features = 10`, which drops barcodes. If
even one barcode was dropped, every score after it landed on the wrong cell —
and the `doublet_scores < 0.4x–0.5x` filters then discarded the wrong cells,
with no warning and no error.

`read_doublet_scores()` in `analysis/helpers.R` closes this two ways:

- a 3-column file (`barcode`, `doublet_score`, `predicted_doublet`) is joined
  **by barcode**, so order and membership stop mattering;
- a legacy 2-column file still works, but only after asserting the row count
  equals the cell count exactly, so a shifted join is a hard error rather than
  quiet corruption.

**Worth checking against your existing data**: run the per-sample analysis with
an old 2-column TSV. If it errors on the row-count assertion, that sample's
published doublet filtering was applied to the wrong cells. If it passes, no
barcodes were dropped and the original result stands.

Regenerate with `scripts/04_scrublet.py` to get the barcoded form. Note the
per-sample `doublet_max` cutoffs in `analysis/config.R` (0.45–0.58) were tuned
against the original Scrublet output, so re-tune them if the new score
distributions differ.

## Known open questions

Two things could not be resolved from the repository and need someone who ran
the original experiments.

### 1. Which reference did the published counts use?

`00_mkref.sh` builds into:

```
.../template_fastas/refdata-gex-GRCh38-2020-A_PV_GFP_mRuby/GRCh38-2020-A_PV_GTF_mRuby
```

which matches the provenance comment in the original analysis script. But the
original count script counted against:

```
.../template_fastas/refdata-gex-GRCh38-2020-A/fasta/GRCh38-2020-A_PV_GTF_mRuby
```

A different parent directory, with a spurious `fasta/` component. At most one
of these is the reference the published matrices came from. `config.sh` uses
the mkref-consistent path; if the other one is correct, change the
`TRANSCRIPTOME` line — but the two would then not be the same object, so
confirm before reusing any existing count matrices.

### 2. Which two samples are the WT GFP replicates?

The original checked for a batch effect by indexing into `Cells()` directly:

```r
Idents(obj, cells[1561:3258])   <- 'Batch1 WT GFP'
Idents(obj, cells[10583:12093]) <- 'Batch2 WT GFP'
```

Those offsets depend on every QC threshold applied above them, so changing any
cutoff would relabel the wrong cells — and the "no significant batch effect"
conclusion would then rest on the wrong comparison. Selection is now by
`orig.ident` via `WT_REPLICATE_SAMPLE_IDS` in `analysis/config.R`.

Given the merge order, `[1561:3258]` (1698 cells) is almost certainly
`WT_GFP_PV`. The second range (1511 cells) falls around the eighth sample,
most plausibly `WT_IRES_GFP_PV`, but that cannot be confirmed without the
per-sample cell counts from that run. **Confirm before relying on the batch
check.**

### 3. The `-2.5` log10 floor is pseudocount-sensitive

`call_infected_status()` fits its mixture on cells whose `log10(percent.virus +
pseudocount)` exceeds `-2.5`. Because the pseudocount is derived from the
smallest non-zero value in each sample, whether zero-virus cells clear that
floor depends on the sample rather than on anything you meant to select for.
This was true of the original too. The function now warns when more than half
the cells entering a fit have zero viral reads, which is the signal that the
floor is not doing what it looks like it does for that sample.

## Tests

```bash
bash tests/run_tests.sh
```

Covers the Slurm memory derivation, shell/Python/R syntax, every R Markdown
chunk, and the two correctness fixes — including a regression test that
reproduces the dropped-cell scenario and asserts the barcode join gives the
right answer where the positional join gave the wrong one.

The R tests need `Seurat` and `mixtools` loadable. On a machine with no package
access, `bash tests/run_tests.sh --stubs` installs the minimal stand-ins in
`tests/stubs/` into a throwaway library. Those stubs reproduce only the two
behaviours under test — `AddMetaData`'s named-vs-unnamed vector contract, and
`normalmixEM` returning components in arbitrary order. Prefer the real packages
where available.

## Data locations

Outputs live on the `ai-fas5.niaid.nih.gov` server under
`lvd_qveu/Lab_Alumni_Archives/Projects_CAM/Scissors`, and in the `CM_kb`
project folder on Skyline.

## Samples

Poliovirus, 14 samples — reporter and polymerase mutants, cotransfections, and
a mock control:

`RFP_C109S_PV`, `WT_GFP_PV`, `WT_GFP_RFP_Y88P_PV`, `WT_GFP_RFP_D177A_PV`,
`RFP_Y88P_PV`, `RFP_D177A_PV`, `WT_GFP_RFP_C109S_PV`, `WT_IRES_GFP_PV`,
`Del_IRES_mRuby3_PV`, `WT_IRES_GFP_Del_IRES_mRuby3_PV`,
`WT_IRES_mRuby3_MutPol_PV`, `Del_IRES_mRuby3_MutPol_PV`,
`WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV`, `Mock_5h_PV`

Plus `CVB3_TT` (coxsackievirus B3) and `EVA71_TT` (enterovirus A71), each
analyzed on its own.
