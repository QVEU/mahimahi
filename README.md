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
| `workflow/Snakefile` | Strand-specific counting workflow (FASTQ or SAM/BAM in). |
| `workflow/scripts/` | `embed_barcodes`, `extract_reads`, `tabulate_strands`, `strand_qc`, `merge_*`. |
| `workflow/envs/scissors.yaml` | Conda environment for the workflow. |
| `config/config.yaml` | Workflow configuration: reference, barcode layout, strand convention, filters. |
| `config/samples.tsv` | One row per sample; input may be FASTQ or SAM/BAM/CRAM. |
| `tests/` | Offline regression tests for the correctness fixes. |
| `tests/workflow/` | End-to-end workflow tests against synthetic ground truth. |

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

## The strand-specific counting workflow

This is the core of SCISSORS: for each cell, count how many positive- and
negative-strand UMIs are logged against each template in the reference. For a
(+)-strand RNA virus the negative strand is the replication intermediate, so
the (-)/(+) ratio measures replication rather than viral load.

It is a Snakemake workflow under `workflow/`, and it takes **either FASTQ or
an already-aligned SAM/BAM/CRAM**, mixed freely in one run.

```
                        config/samples.tsv
                               |
        +----------------------+----------------------+
        |                      |                      |
   FASTQ, barcode        FASTQ, barcode          SAM/BAM/CRAM
   in the read name      in R1 (raw 10x)         (already aligned)
        |                      |                      |
        |              embed_barcodes                 |
        |              (R1 -> R2 name)                |
        |                      |                      |
        +------- align --------+            normalise_alignments
             (minimap2 -ax sr)                        |
                    |                                 |
                    +----------------+----------------+
                                     |
                              extract_reads
                    per read: CBC, UMI, ref_name, strand
                                     |
                    +----------------+----------------+
                    |                                 |
            tabulate_strands                     strand_qc
   distinct UMIs per cell x template      strand balance per template,
        x strand; Rep_Index                asserts the convention
                    |                                 |
              merge_counts                    merge_strand_qc
        results/scissors_counts.tsv.gz   results/strand_qc_summary.tsv
```

### Running it

```bash
# describe your samples and reference
$EDITOR config/samples.tsv config/config.yaml

snakemake --cores 8                    # run
snakemake --cores 8 -n                 # dry run: show the plan
snakemake --cores 8 --software-deployment-method conda   # managed deps
```

`config/samples.tsv` needs `sample` and `input` per row. The input's extension
decides the route: a FASTQ gets aligned against `template`, a SAM/BAM/CRAM is
used as-is. Optional columns are `mate_fastq` (the R1 barcode read),
`barcode_source`, `datalabel`, and `whitelist`.

### Where the barcode comes from

Set `barcode.source` in `config/config.yaml`, or per sample:

| source | use for | barcode read from |
| --- | --- | --- |
| `read_name` | DRAGEN scRNA output | a colon-delimited field of the read name |
| `tags` | Cell Ranger BAM | `CB`/`UB` tags |
| `paired_fastq` | raw 10x FASTQ pair | first bases of R1, embedded into the R2 name |

`barcode_length` and `umi_length` are configuration, not literals: **10x 3' v2
is 16+10, v3 and v3.1 are 16+12.** `mahimahi_dragen.py` hardcoded a 10 nt UMI,
which on v3 discards two bases and shrinks the UMI space 16-fold. The
resulting undercount is under 1% at typical depth, but it lands mostly on
`Pos` and so inflates `Neg/Pos` systematically. Set these to the kit you ran.

### The strand convention, and how to check it

Mapped orientation is not the strand of the original RNA -- the relationship
depends on the chemistry and on which mate was aligned. Inverting it swaps
`Pos` and `Neg` on every template and produces output that looks perfectly
reasonable.

`strand.convention` defaults to `reverse_is_positive`, which is what every
existing SCISSORS script used and what the published output corroborates
(`Pos` dominant, `Neg` at 0.7-4%).

Rather than trust it, set `strand.sense_control` to a reference whose
orientation you know -- a host transcript, a spike-in, a reporter cassette.
`rule strand_qc` then requires that sequence to come out predominantly
positive-sense and **fails the run** if it does not:

```
strand_qc FAILED:
  - sense control 'HostControl' is only 0.0% positive-sense (require >=90%).
    The strand convention is probably inverted.
```

### What this fixes

The workflow replaces `mahimahi.sh`, `mahimahi_dragen.py`, and the four
divergent copies of the `SCISSORS()` tabulation. The substantive changes:

- **UMIs are deduplicated within (CBC, ref_name, strand)**, and the wide table
  is keyed on `(CBC, ref_name)`. Both Python versions of `SCISSORS()` used
  `pivot_table(index="CBC")`, whose default `aggfunc` is `"mean"` -- so `Neg`,
  `Pos`, `Rep_Index` and `Neg_PosRatio` were averaged across every template in
  a cell, weighted by read count. `Rep_Index` came out identical for every
  template, so it could not distinguish donor from acceptor replication. See
  `NEWS.md` for the worked example.
- **Alignments are selected by FLAG bits**, not `flag in (0, 16)`. The original
  test is right for single-end minimap2 output but discards nearly everything
  in a paired or Cell Ranger BAM.
- **Read-name layout is validated up front** against a sample of reads, with an
  error naming the offending read, instead of raising `IndexError` partway
  through a file.
- **Output columns are selected by name**, not `.iloc[:, 12:]` and `.iloc[:, 3]`.
- **Output goes where it is told.** The original wrote to the interpreter's
  working directory, so two input directories each holding a
  `Mock_5h_S5_CBC.csv` silently overwrote one another.
- **One filter, in one place.** The originals filtered `CBC_readcount > 100` in
  Python and then `UMI_count > 100` in R -- two different quantities behind the
  same number.
- **R1 is no longer aligned pointlessly.** `mahimahi.sh` ran minimap2 on the
  barcode read against the viral template, and the script consuming the output
  ignored that SAM entirely.

### Testing it

```bash
bash tests/workflow/run_workflow_tests.sh
```

Generates synthetic reads whose per-cell, per-template, per-strand UMI counts
are known exactly, runs all three input routes, and checks every value. The
fixtures deliberately give one cell two templates with different replication
levels (eGFP 0.0909, mRuby3 0.3333), which is the case the original collapsed;
the suite asserts they come out distinct. It also confirms that an inverted
strand convention is rejected, and that a second run is a no-op.


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
