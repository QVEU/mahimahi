# Legacy analyses — retained for provenance, not for use

These notebooks are the original SCISSORS analyses, kept unchanged because they
are the record of what produced the results so far. **Do not run them.**
Their replacements are listed below, and [../docs/MIGRATION.md](../docs/MIGRATION.md)
explains each change and why.

They are archived here because they exist nowhere else. The four original
mahimahi scripts are *not* kept: they lived in the original mahimahi
repository (commit `4c223b9`), which was deleted in October 2026 when this
repository took over the mahimahi name.

## What is here

| Notebook | Language | What it did | Replaced by |
| --- | --- | --- | --- |
| `process_mahimahi.ipynb` | Python 3.9 | `SCISSORS()`: per-read CSV → per-cell strand counts. Grouped by `("CBC","flag")`, **omitting `ref_name`**, so no per-template resolution at all. | `workflow/scripts/tabulate_strands.py` |
| `process_mahimahi-Working020724.ipynb` | Python 3.9 | Same, with `ref_name` added to the groupby — a partial fix. The `pivot_table(index="CBC")` was left alone, so `Neg`/`Pos`/`Rep_Index` remained averaged across templates. | `workflow/scripts/tabulate_strands.py` |
| `mahimahi.Rnb.ipynb` | R | The **only correct tabulation**: `dcast(CBC+strand ~ ref_name, fun.aggregate = length(unique(UMI)))`, plus DRAGEN barcode-whitelist filtering. Its logic is what the workflow implements. | `workflow/scripts/tabulate_strands.py`, whitelist via the `whitelist` column in `config/samples.tsv` |
| `mahimahi_81423.ipynb` | R 4.3.1 | `allreadsdir()` does `spread(UMI, UMI)`, which would make one column per UMI *value*. Abandoned. | — |
| `mahimahi_Analysis_Plotting.ipynb` | R | Earlier plotting pass over the strand counts. | `analysis/mahimahi_replication.R` |
| `Scissors_Analysis_v4.ipynb` | R (mislabelled `Python 3`) | The main replication analysis and figures. | `analysis/mahimahi_replication.R` + `analysis/replication.R` |

Cell outputs are left embedded. They are most of the file size, and they are
also the evidence of what these notebooks actually produced.

## Why they should not be run

Each of these has at least one defect that a reader should know about before
treating its output as correct.

**`Scissors_Analysis_v4.ipynb` — `fitSet()` could not return.** It has a bare
`break` between its two `Fits` assignments. In R, `break` outside a loop is an
error (`no loop for break/next, jumping to top level`), so the call never
reached the block that built the `Ratio`/`Error` frame. Every figure depending
on those columns — cells 9, 16, 19, 20, 21, 24 — rested on a call that could
not complete. It also extracts coefficients positionally
(`summary(glm(...))$coefficients[2]` and `[4]`), which is only correct while
the coefficient matrix keeps its shape, and plots `Ratio ± Error` where `Error`
is the slope's standard error — a ~68% interval drawn as though it were 95%.

**Both Python `SCISSORS()` versions averaged across templates.**
`pivot_table(index="CBC", ...)` defaults `aggfunc` to `"mean"`, and
`UMI_strand_count` is repeated once per read, so `Neg` and `Pos` came out as
read-count-weighted averages over every template in the cell. Worked example,
one cell:

```
truth     eGFP   Pos=1000 Neg=50   Rep_Index=0.0476
          mRuby3 Pos= 100 Neg=20   Rep_Index=0.1667
produced  both templates  Neg=41.43  Pos=918.18  Rep_Index=0.0432
```

`Rep_Index` therefore could not distinguish donor from acceptor replication.

**The notebook kernels are mislabelled.** `Scissors_Analysis_v4.ipynb`
declares a `Python 3` kernel but contains R throughout; opening it with the
declared kernel fails on the first cell.

## Which produced the existing figures?

This has not been established, and the three tabulations do not agree with each
other. None of it is published, so the fix is to regenerate through `workflow/`
rather than to work out retrospectively which notebook was used; the test below
is only worth running if you need to know whether an existing figure was wrong.

`Scissors_Analysis_v4.ipynb` cell 19 plots per-cell `Rep_Index` for eGFP
against mRuby3. Under the `pivot_table` bug both templates carry the same
number, so **that plot must be a perfect y = x line**. If it is, the counts
came from a Python version and the donor/acceptor result needs regenerating
through `workflow/`. If the points scatter, they came from
`mahimahi.Rnb.ipynb` and stand.

`analysis/mahimahi_replication.R` checks this automatically and warns when more
than 99% of co-infected cells share an identical `Rep_Index` across templates.

## Package versions these assumed

Recorded because none of them are current, and two are no longer installable
as specified:

- **Seurat 4.4.0 / SeuratObject 4.1.4** — both archived on CRAN. The
  replacements target Seurat 5, where `merge()` is layer-aware.
- **standalone `scrublet` 0.2.3** — unmaintained; superseded by
  `scanpy.pp.scrublet`.
- **Cell Ranger 7.2.0**, **bcl2fastq2 2.20** — bcl2fastq2 is end-of-life at
  Illumina, and `--create-bam` became mandatory for `cellranger count` in 8.0.
- Python 3.9 for the `process_mahimahi` notebooks; R 4.3.1 for
  `mahimahi_81423.ipynb`.
