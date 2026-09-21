# Changes from QVEU/SCISSORS

This repository restructures the original pipeline. Analysis intent is
unchanged; the changes below fix defects, remove hardcoded paths and indices,
and make the shell scripts runnable as submitted. Hand-tuned QC thresholds
were carried over verbatim.

## Correctness

These change results, or would have under some conditions.

**Doublet scores could be attached to the wrong cells.**
`AddMetaData` given an unnamed vector assigns by position. Scores came from
2-column TSVs with no barcodes, computed on the full
`filtered_feature_bc_matrix`, while the Seurat objects were built with
`min.features = 10`, which drops barcodes. One dropped barcode shifted every
subsequent score onto the wrong cell, and the `doublet_scores` filters then
discarded the wrong cells silently. `read_doublet_scores()` joins by barcode
where available and hard-errors on a length mismatch otherwise.
Affected all three analyses. See README "Doublet scores".

**Infected/uninfected labels depended on mixtools' component order.**
The loop assumed `normalmixEM` returned `mu[1]` as the uninfected mode and
`mu[2]` as the high-replication mode. It makes no such guarantee — component
order follows random initialization. A run returning them swapped would invert
every call in that sample with nothing visibly wrong. Components are now
sorted by `mu`, and the seed is fixed. Affected
`SCISSORS_SeuratAnalysis.r`.

**The mixture was fit on one scale and selected on another.**
Cells were selected with a `log10` cutoff but the mixture was fit on raw
`percent.virus`. The fit is now on `log10` throughout, with thresholds
converted back only for labelling.

**`break` where `next` was meant.**
The infected-calling loop used `break` on the mock sample, which exits the loop
rather than skipping one iteration. Harmless only because `Mock_5h_PV` was last
in the vector; reordering it would have left every later sample without an
`InfectedStatus`, with no error.

**Batch assignment by hardcoded cell indices.**
`cells[1561:3258]` and `cells[10583:12093]` labelled the WT GFP replicates.
Those offsets are a function of every QC threshold above them, so any cutoff
change would relabel the wrong cells and the "no significant batch effect"
conclusion would rest on the wrong comparison. Selection is now by
`orig.ident`. See README "Known open questions" — the second sample needs
confirming.

**Volcano plot painted over its own data.**
Grey "not significant" points were drawn *after* the coloured UP/DOWN layers,
hiding them. Grey is now the base layer. The y axis was also clipped at 30000
via `coord_cartesian`; that clip is removed.

**Degenerate mixture inputs now error.**
A fit set with fewer than 5 distinct values converges to two identical
zero-variance components and calls everything above a single point infected.
This is now a hard error. A warning also fires when more than half the cells
entering a fit have zero viral reads, which indicates the `-2.5` log floor is
being swamped by the pseudocount.

**`half`min was not half the minimum.**
`halfmin_virus` was `min(positive)`, despite the name. Now `min(positive) / 2`,
the conventional pseudocount.

## Scripts that could not run as submitted

**`cellranger count.sh` had a broken shebang.** `#!/bin/bash#` — the trailing
`#` is part of the interpreter path, so direct execution failed.

**Two scripts overcommitted memory and would be OOM-killed.**
`mkfastq` paired `#SBATCH --mem=12G` with `--localmem=160`; `aggregate` paired
`--mem=10G` with `--localmem 100`. Cell Ranger schedules work believing it has
memory the cgroup will not grant. `--localcores`/`--localmem` are now derived
from the Slurm allocation in `config.sh`, so the mismatch cannot recur.

**Invalid Slurm directives.** `#SBATCH - cwd` is SGE's `-cwd`;
`#SBATCH --ntasks-per-node = 16` has spaces around the `=`, which Slurm
rejects. Both removed.

**`--ntasks-per-*` did not reserve the requested CPUs.** All scripts used
`--ntasks-per-node` or `--ntasks-per-core` alongside `--localcores=16`. Now
`--cpus-per-task=16`.

**`MkRefs.sh` blocked on an interactive shell.** `srun --pty bash` on line 3
meant the `mkref` call below it did not run until you exited that shell, and
then ran on the login node rather than in the allocation. Now a plain batch
script.

**`count.sh` could only ever process one sample.** `--sample`, `--id`,
`--fastqs` and `--output-dir` were hardcoded to `Mock_5h`, and it accepted a
`runID` argument it never used — so 13 of 14 samples were not reproducible from
the repository. Parameters now come from `scripts/samples.tsv`, with
Slurm-array support.

**Added preflight checks.** Missing transcriptome, missing FASTQ directory, no
FASTQs matching the sample prefix, unfilled `TODO` rows, and existing output
directories are all caught before Cell Ranger starts. Counting zero reads
otherwise succeeds and yields an empty matrix noticed only much later.

## Reproducibility

**`.r` files that were R Markdown.** `CVB3.r` and `EVA71.r` had YAML front
matter and knitr chunks — `Rscript CVB3.r` failed on line 1. Now `.Rmd`.

**Paths.** The same share was referred to four ways: `/Volumes/lvd_qve` (10x),
`/Volumes/LVD_QVE` (2x — a case difference that breaks on case-sensitive
volumes), `/data/lvd_qve` (4x), and `~/lab_share` (2x). One script read from
`~/lab_share` and wrote to `/data/lvd_qve` in the same run. `analysis/config.R`
resolves the share once, with a `SCISSORS_SHARE_ROOT` override.

**Package installation on every run.** `remotes::install_version()` ran
unconditionally at the top of the analysis script, reinstalling Seurat on every
`source()`. Now a documented one-time step, with a Seurat v4 version assertion.

**Dead code.** `nestorawa_forcellcycle_expressionMatrix.txt` was read into
`exp.mat` in all three analyses and never used — the cell cycle genes come from
Seurat's built-in `cc.genes`. Two of the three read it from `~/Downloads/`.
Removed.

**85 lines of copy-pasted doublet blocks** became one loop.

**Doublet filename disagreement.** The per-sample analyses read
`<SAMPLE>_Doublet_scores.tsv`; the integrated one read a bare
`Doublet_scores.tsv`. Both are now accepted.

**`ncol` mismatches in `VlnPlot`.** `ncol = 5` was passed for 3 and 4 features.

**Missing scrublet stage.** `scrublet.Rmd` was referenced in comments but never
committed, so the doublet scores could not be regenerated.
`scripts/04_scrublet.py` is a reimplementation using Scrublet's documented
defaults — **its parameters are not recovered from the original**, so validate
against existing TSVs before reusing it on published data.

**Misleading section header.** "Perform integration analysis" preceded a plain
`merge()` with no integration or batch correction. Retitled.

**Missing input documentation.** `SampleSheet.csv` and the `aggr` CSV were
undocumented inputs. Their formats are now in the README.

## Added

- `config.sh` / `analysis/config.R` — paths, sample lists, QC thresholds.
- `analysis/helpers.R` — doublet joining, infected calling, QC filtering.
- `scripts/samples.tsv` — per-sample FASTQ paths.
- `tests/` — offline regression tests, including one reproducing the
  dropped-cell doublet scenario.
- `.gitignore`.
