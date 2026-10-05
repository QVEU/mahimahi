# Where the old scripts went

SCISSORS accumulated several implementations of the same steps, which
disagreed with each other. This maps every original file to its replacement
and says why.

**Where the originals are.** The four mahimahi scripts live in
[QVEU/mahimahi-legacy](https://github.com/QVEU/mahimahi-legacy) at commit `4c223b9`, verified
byte-identical to the versions reviewed here; that repository is deprecated and
archived. The six notebooks had no home anywhere and are kept in
[`../legacy/`](../legacy/README.md), unchanged, with their outputs embedded.

## Retired

| Original | Replacement | Why |
| --- | --- | --- |
| `mahimahi.sh` | `workflow/Snakefile` (`rule align`) | `#!/bin/bash/` (trailing slash) is not a valid interpreter path. Aligned R1 against the viral template and then called a script that ignored the R1 SAM. |
| `mahimahi_dragen.sh` | `workflow/Snakefile` (`rule align`) | Superseded; the workflow handles FASTQ and pre-aligned input. |
| `mahimahi.py` | `workflow/scripts/extract_reads.py` | **Never ran.** Calls `stop()`, which does not exist in Python, and `pd.DataFrame(readDict[i])` on a flat dict of scalars raises `ValueError`. |
| `mahimahi_dragen.py` | `workflow/scripts/extract_reads.py` | Hardcoded `[0:16]`/`[16:26]` (truncates a v3 UMI); one row per read name only; `flag in (0,16)`; wrote seq/qual to CSV; overwrote the input SAM when the filename did not match the expected pattern. |
| `process_mahimahi.ipynb` | `workflow/scripts/tabulate_strands.py` | `groupby(["CBC","flag"])` omits `ref_name`, so no per-template resolution at all. |
| `process_mahimahi-Working020724.ipynb` | `workflow/scripts/tabulate_strands.py` | Adds `ref_name` to the groupby, but still `pivot_table(index="CBC")`, so `Neg`/`Pos`/`Rep_Index`/`Neg_PosRatio` are read-count-weighted averages across templates. |
| `mahimahi.Rnb.ipynb` | `workflow/scripts/tabulate_strands.py` | The only *correct* tabulation: `dcast(CBC+strand ~ ref_name, fun.aggregate = length(unique(UMI)))`, plus DRAGEN whitelist filtering. Its logic is what the workflow implements; the whitelist filter is the `whitelist` column in `config/samples.tsv`. |
| `mahimahi_81423.ipynb` | — | `allreadsdir()` does `spread(UMI, UMI)`, which would make one column per UMI value. Abandoned. |
| `Scissors_Analysis_v4.ipynb` | `analysis/scissors_replication.R` + `analysis/replication.R` | `fitSet()` could not return (bare `break`); `CollectFiles()` dropped the `Neg`/`Pos` columns a later cell needed; stale factor levels silently dropped a sample; axis labels described `Rep_Index` while plotting the slope; several cells referenced undefined objects. |
| `SCISSORS_SeuratAnalysis.r` | `analysis/PV_mutants_integrated.R` | See `NEWS.md`. |
| `CVB3.r`, `EVA71.r` | `analysis/CVB3_QC.Rmd`, `analysis/EVA71_QC.Rmd` | Were R Markdown with a `.r` extension, so `Rscript` could never run them. |

## Which version produced your existing CSVs?

Nothing here has been published, so this is a question about what your current
working numbers mean, not about whether a result has to be corrected. The
cheapest answer is usually to skip the forensics: regenerate through
`workflow/`, which is tested against synthetic ground truth, and compare. The
diagnostic below is worth running only if regenerating is expensive -- the raw
SAM/BAMs are gone, say, or a figure is already in a talk and you want to know
whether it was wrong.

It matters because the two Python tabulations and the R one do not agree.

`Scissors_Analysis_v4.ipynb` cell 19 plots per-cell `Rep_Index` for eGFP
against mRuby3:

```r
ggplot(dcast.data.table(PVList, CBC+sample+filename ~ ref_name,
                        value.var = "Rep_Index", fun.aggregate = max)) +
    geom_point(aes(eGFP, mRuby3))
```

Under the `pivot_table(index="CBC")` bug both templates carry the same number,
so **that plot must be a perfect y=x line**. If it is, the counts came from a
Python version and the donor/acceptor result needs regenerating through the
workflow. If the points scatter, they came from `mahimahi.Rnb.ipynb` and stand.

`analysis/scissors_replication.R` checks this automatically and warns when more
than 99% of co-infected cells have identical `Rep_Index` across templates.

## Column mapping

The workflow's `results/scissors_counts.tsv.gz` keeps the original column
names, so existing analysis code mostly ports directly.

| `*_CBC_SCISSORS.csv` | `scissors_counts.tsv.gz` | Note |
| --- | --- | --- |
| `CBC` | `CBC` | same |
| `ref_name` | `ref_name` | same |
| `sample` | `sample` | was the mahimahi *filename*; now the sample name from the sheet |
| `datalabel` | `datalabel` | same derivation (`_S<n>` stripped), but only a genuine `_S<digits>` tail |
| `CBC_readcount` | `CBC_readcount` | same |
| `UMI_count` | `UMI_count` | distinct viral UMIs per cell |
| `Neg`, `Pos` | `Neg`, `Pos` | **now per (cell, template)**, previously averaged across templates |
| `Neg_PosRatio` | `Neg_PosRatio` | same formula, `Neg/(Pos+1)`; now per template |
| `Rep_Index` | `Rep_Index` | same formula, `Neg/(Pos+Neg)`; now per template |
| `det_strand`, `UMI_strand_count` | in `<sample>/counts_long.tsv.gz` | the long form is kept per sample |

## `Ratio` is not `Rep_Index`

`fitSet()` returned `Ratio`, the slope of `Neg ~ Pos` across cells. That is
(−)/(+). `Rep_Index` is (−)/total. The notebook's axes said
`((-)strand/total(vRNA))` over plots of `Ratio`.

They converge only as the ratio approaches zero, since (−)/total = r/(1+r):

| slope r = (−)/(+) | (−)/total | difference |
| --- | --- | --- |
| 0.01 | 0.0099 | 1% |
| 0.05 | 0.0476 | 5% |
| 0.20 | 0.1667 | 17% |

`analysis/replication.R` exposes both and labels each for what it is:
`fit_replication_slope()` returns `slope` with a 95% CI, and
`summarise_rep_index()` returns `mean_rep_index`/`median_rep_index`.
The originals plotted `Ratio ± Error` where `Error` was the standard error of
the slope — roughly a 68% interval drawn as though it were 95%.
