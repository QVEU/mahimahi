## Replication-rate estimation from strand-specific counts.
##
## Pure computation: no plotting, no data.table, no tidyverse, so it can be
## unit-tested and used from any script. Figures live in
## analysis/scissors_replication.R.
##
## Replaces fitSet() from Scissors_Analysis_v4.ipynb.

## ---------------------------------------------------------------------------
## Two different quantities, kept apart
##
## The notebook used one label for both, which is where the axis confusion came
## from. They answer different questions and are not interchangeable:
##
##   rep_index  Per cell, per template: Neg / (Pos + Neg). The fraction of that
##              template's UMIs in a given cell that are negative-sense. This is
##              literally "(-)strand / total(vRNA)" -- the label the notebook
##              put on its axes. Computed by the workflow.
##
##   slope      Per group (sample x template): the slope of a linear fit of Neg
##              on Pos across cells. This is (-)/(+), NOT (-)/total. It is a
##              population-level estimate that is robust to per-cell depth,
##              which is why the original used it. But an axis reading
##              "(-)strand/total(vRNA)" over a slope is mislabelled; at small
##              ratios the two are close, but they are not the same number.
##
## For small r, (-)/total = r/(1+r) where r = (-)/(+), so slope and mean
## rep_index converge as r -> 0 and diverge visibly above r ~ 0.1.
## ---------------------------------------------------------------------------

REQUIRED_COUNT_COLUMNS <- c("CBC", "ref_name", "sample", "Pos", "Neg")


check_counts <- function(counts) {
  missing <- setdiff(REQUIRED_COUNT_COLUMNS, colnames(counts))
  if (length(missing) > 0) {
    stop("strand count table is missing required column(s): ",
         paste(missing, collapse = ", "),
         "\nExpected the output of the SCISSORS workflow ",
         "(results/scissors_counts.tsv.gz).")
  }
  invisible(TRUE)
}


## ---------------------------------------------------------------------------
## fit_replication_slope
##
## Fixes relative to fitSet():
##
## 1. It returns. fitSet() had a bare `break` between its two Fits assignments.
##    `break` outside a loop is an error in R -- "no loop for break/next,
##    jumping to top level" -- so every call failed before reaching the block
##    that built the Ratio/Error data frame. Verified against R 4.3.3. Every
##    figure downstream of it (Ratio, Error) was drawing on a call that could
##    not complete.
##
## 2. Coefficients are extracted by name. The original used
##    summary(glm(...))$coefficients[2] and [4], relying on column-major
##    flattening of a 2x4 matrix. Correct for a well-behaved fit, but if a
##    group has too few cells or no variance in Pos the matrix changes shape
##    and those indices silently return the wrong number.
##
## 3. Groups that cannot be fitted are reported, not fatal. One sample with
##    three cells no longer takes down the whole run.
##
## 4. A confidence interval, not a bare standard error. The original plotted
##    Ratio +/- Error where Error was the standard error of the slope, i.e.
##    roughly a 68% interval drawn as though it were a 95% one.
##
## `min_umis` replaces the hardcoded UMI_count > 100 in the notebook and the
## CBC_readcount > 100 in the Python. Those were two different quantities
## behind the same number, applied in sequence; filter once, explicitly.
## ---------------------------------------------------------------------------
fit_replication_slope <- function(counts,
                                  by = c("sample", "ref_name"),
                                  min_umis = 0,
                                  min_cells = 5,
                                  conf_level = 0.95,
                                  exclude_refs = character()) {
  check_counts(counts)
  missing_by <- setdiff(by, colnames(counts))
  if (length(missing_by) > 0) {
    stop("grouping column(s) not in the table: ", paste(missing_by, collapse = ", "))
  }

  if (length(exclude_refs) > 0) {
    counts <- counts[!counts$ref_name %in% exclude_refs, , drop = FALSE]
  }

  if (min_umis > 0) {
    if (!"UMI_count" %in% colnames(counts)) {
      stop("min_umis was requested but the table has no UMI_count column.")
    }
    counts <- counts[!is.na(counts$UMI_count) & counts$UMI_count >= min_umis, ,
                     drop = FALSE]
  }

  if (nrow(counts) == 0) {
    warning("no rows left after filtering; returning an empty result.",
            call. = FALSE)
    return(empty_slope_frame(by))
  }

  groups <- interaction(counts[by], drop = TRUE, sep = "\r")
  pieces <- split(seq_len(nrow(counts)), groups)

  rows <- lapply(names(pieces), function(key) {
    idx <- pieces[[key]]
    labels <- as.list(strsplit(key, "\r", fixed = TRUE)[[1]])
    names(labels) <- by
    fit_one_group(counts[idx, , drop = FALSE], labels, min_cells, conf_level)
  })

  result <- do.call(rbind, c(rows, list(make.row.names = FALSE)))
  rownames(result) <- NULL

  failed <- result$status != "ok"
  if (any(failed)) {
    message(sum(failed), " of ", nrow(result),
            " groups could not be fitted (see the status column): ",
            paste(unique(result$status[failed]), collapse = ", "))
  }
  result
}


empty_slope_frame <- function(by) {
  base <- lapply(by, function(x) character())
  names(base) <- by
  do.call(data.frame, c(base, list(
    n_cells = integer(), slope = numeric(), slope_se = numeric(),
    conf_low = numeric(), conf_high = numeric(), intercept = numeric(),
    r_squared = numeric(), status = character(),
    stringsAsFactors = FALSE)))
}


fit_one_group <- function(group, labels, min_cells, conf_level) {
  na_row <- function(status, n_cells) {
    do.call(data.frame, c(labels, list(
      n_cells = n_cells, slope = NA_real_, slope_se = NA_real_,
      conf_low = NA_real_, conf_high = NA_real_, intercept = NA_real_,
      r_squared = NA_real_, status = status,
      stringsAsFactors = FALSE)))
  }

  usable <- !is.na(group$Pos) & !is.na(group$Neg)
  group <- group[usable, , drop = FALSE]
  n <- nrow(group)

  if (n < min_cells)                    return(na_row("insufficient_cells", n))
  if (length(unique(group$Pos)) < 2)    return(na_row("no_variance_in_Pos", n))
  if (all(group$Neg == 0))              return(na_row("no_negative_strand", n))

  fit <- try(stats::lm(Neg ~ Pos, data = group), silent = TRUE)
  if (inherits(fit, "try-error"))       return(na_row("fit_failed", n))

  coefs <- stats::coef(summary(fit))
  if (!"Pos" %in% rownames(coefs))      return(na_row("no_slope_term", n))

  ## Named indexing, so a changed matrix shape cannot quietly hand back the
  ## intercept's standard error in place of the slope's.
  slope <- coefs["Pos", "Estimate"]
  slope_se <- coefs["Pos", "Std. Error"]
  intercept <- if ("(Intercept)" %in% rownames(coefs))
    coefs["(Intercept)", "Estimate"] else NA_real_

  df_resid <- stats::df.residual(fit)
  half_width <- if (is.finite(df_resid) && df_resid > 0) {
    stats::qt(1 - (1 - conf_level) / 2, df = df_resid) * slope_se
  } else NA_real_

  do.call(data.frame, c(labels, list(
    n_cells = n,
    slope = slope,
    slope_se = slope_se,
    conf_low = slope - half_width,
    conf_high = slope + half_width,
    intercept = intercept,
    r_squared = summary(fit)$r.squared,
    status = "ok",
    stringsAsFactors = FALSE)))
}


## ---------------------------------------------------------------------------
## Per-cell replication index summaries.
##
## rep_index is computed per (cell, template) by the workflow. This summarises
## it without going through a regression, which is the right thing to plot on
## an axis labelled "(-)strand / total(vRNA)".
## ---------------------------------------------------------------------------
summarise_rep_index <- function(counts,
                                by = c("sample", "ref_name"),
                                min_umis = 0,
                                exclude_refs = character()) {
  check_counts(counts)
  if (!"Rep_Index" %in% colnames(counts)) {
    stop("table has no Rep_Index column; expected the SCISSORS workflow output.")
  }
  if (length(exclude_refs) > 0) {
    counts <- counts[!counts$ref_name %in% exclude_refs, , drop = FALSE]
  }
  if (min_umis > 0) {
    counts <- counts[!is.na(counts$UMI_count) & counts$UMI_count >= min_umis, ,
                     drop = FALSE]
  }

  groups <- interaction(counts[by], drop = TRUE, sep = "\r")
  pieces <- split(seq_len(nrow(counts)), groups)

  rows <- lapply(names(pieces), function(key) {
    idx <- pieces[[key]]
    labels <- as.list(strsplit(key, "\r", fixed = TRUE)[[1]])
    names(labels) <- by
    values <- counts$Rep_Index[idx]
    values <- values[!is.na(values)]
    do.call(data.frame, c(labels, list(
      n_cells = length(values),
      mean_rep_index = if (length(values)) mean(values) else NA_real_,
      median_rep_index = if (length(values)) stats::median(values) else NA_real_,
      q25 = if (length(values)) unname(stats::quantile(values, 0.25)) else NA_real_,
      q75 = if (length(values)) unname(stats::quantile(values, 0.75)) else NA_real_,
      stringsAsFactors = FALSE)))
  })

  result <- do.call(rbind, c(rows, list(make.row.names = FALSE)))
  rownames(result) <- NULL
  result
}


## ---------------------------------------------------------------------------
## Joining strand counts to Seurat metadata
##
## The notebook did:
##   CVBsc[, CBC := tstrsplit(V1, "-")[1]]
##   merge.data.table(scData, InputSet, by = c("CBC", "datalabel"))
##
## an inner join on a bare string. Two things make that lossy in ways nothing
## reported:
##
##   - The workflow's barcodes come from the read as sequenced. Cell Ranger's
##     are error-corrected against the 10x whitelist. A barcode with a
##     sequencing error is corrected on one side and not the other, so the cell
##     silently fails to join.
##   - An inner join drops non-matching rows from both sides without comment,
##     so a join that loses a third of the cells looks the same as one that
##     loses none.
##
## This reports both directions and can fail when the loss exceeds a
## threshold, so a bad join is visible rather than inferred from a thin plot.
## ---------------------------------------------------------------------------
join_seurat_metadata <- function(counts, metadata,
                                 counts_key = "CBC",
                                 metadata_key = "cell_barcode",
                                 by_extra = character(),
                                 max_unmatched_frac = 0.25,
                                 on_excess = c("warn", "error")) {
  on_excess <- match.arg(on_excess)
  for (nm in c(counts_key, by_extra)) {
    if (!nm %in% colnames(counts)) stop("counts has no column '", nm, "'")
  }
  for (nm in c(metadata_key, by_extra)) {
    if (!nm %in% colnames(metadata)) stop("metadata has no column '", nm, "'")
  }

  counts$.join_key <- strip_barcode_suffix(counts[[counts_key]])
  metadata$.join_key <- strip_barcode_suffix(metadata[[metadata_key]])

  keys <- c(".join_key", by_extra)
  merged <- merge(counts, metadata, by = keys, all = FALSE, sort = FALSE)

  counts_cells <- unique(counts$.join_key)
  meta_cells <- unique(metadata$.join_key)
  matched <- unique(merged$.join_key)

  unmatched_counts <- setdiff(counts_cells, matched)
  unmatched_meta <- setdiff(meta_cells, matched)
  frac <- if (length(counts_cells)) length(unmatched_counts) / length(counts_cells) else 0

  message(sprintf(
    paste0("join: %d/%d strand-count cells matched Seurat metadata (%.1f%% unmatched); ",
           "%d of %d Seurat cells had no strand counts"),
    length(matched), length(counts_cells), 100 * frac,
    length(unmatched_meta), length(meta_cells)))

  if (frac > max_unmatched_frac) {
    msg <- sprintf(
      paste0("%.1f%% of strand-count cells did not match any Seurat cell ",
             "(threshold %.0f%%).\n",
             "  Likely causes: barcodes error-corrected on one side only, a ",
             "barcode suffix mismatch, or mismatched samples.\n",
             "  Example unmatched barcodes: %s"),
      100 * frac, 100 * max_unmatched_frac,
      paste(utils::head(unmatched_counts, 5), collapse = ", "))
    if (on_excess == "error") stop(msg) else warning(msg, call. = FALSE)
  }

  merged$.join_key <- NULL
  attr(merged, "join_report") <- list(
    matched_cells = length(matched),
    counts_cells = length(counts_cells),
    metadata_cells = length(meta_cells),
    unmatched_counts = unmatched_counts,
    unmatched_metadata = unmatched_meta,
    unmatched_fraction = frac)
  merged
}


## Drop a 10x '-1' or merged '-<sample>' suffix.
strip_barcode_suffix <- function(x) sub("-.*$", "", as.character(x))


## ---------------------------------------------------------------------------
## read_counts_table: read a possibly-gzipped TSV without requiring R.utils
##
## data.table::fread() below ~1.15 refuses .gz unless R.utils is installed,
## which it often is not on a cluster R. Streaming through gzip keeps memory
## bounded; the readLines path is a fallback for machines without it.
## ---------------------------------------------------------------------------
read_counts_table <- function(path) {
  if (!file.exists(path)) stop("file not found: ", path)
  gzipped <- grepl("\\.gz$", path)

  if (requireNamespace("data.table", quietly = TRUE)) {
    fread <- data.table::fread
    if (!gzipped) return(as.data.frame(fread(path)))
    if (nzchar(Sys.which("gzip"))) {
      return(as.data.frame(fread(cmd = paste("gzip -dc", shQuote(path)))))
    }
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    return(as.data.frame(fread(text = readLines(con))))
  }

  ## No data.table at all: base R handles gz connections directly.
  con <- if (gzipped) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, stringsAsFactors = FALSE)
}
