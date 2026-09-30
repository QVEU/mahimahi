## Shared helpers for the SCISSORS R analyses.
## Sourced after analysis/config.R.

## ---------------------------------------------------------------------------
## read_doublet_scores: attach Scrublet scores to a Seurat object by barcode
##
## Why this exists
## ---------------
## The pre-restructure scripts did:
##
##   doublets <- read.table(path, header = FALSE)
##   obj <- AddMetaData(obj, metadata = doublets$doublet_scores, ...)
##
## AddMetaData given an unnamed vector assigns it to cells in object order. The
## Scrublet files were computed on the full filtered_feature_bc_matrix, but the
## Seurat objects were built with min.features = 10, which drops barcodes. If
## even one barcode was dropped, every score after it landed on the wrong cell,
## and the doublet_scores < 0.4x-0.5x filters downstream then discarded the
## wrong cells -- with no warning and no error.
##
## This helper removes the failure mode two ways:
##   - a 3-column file (barcode, score, predicted) is joined by barcode, so
##     order and membership stop mattering entirely;
##   - a legacy 2-column file still works positionally, but only after
##     asserting the row count matches the cell count exactly, so a shifted
##     join becomes a hard error instead of quiet corruption.
##
## Regenerate scores with scripts/04_scrublet.py to get the barcoded form.
## ---------------------------------------------------------------------------
read_doublet_scores <- function(object, path, sample_id = NULL) {
  if (!file.exists(path)) {
    stop("Doublet score file not found: ", path,
         "\nRun scripts/04_scrublet.py for this sample first.")
  }

  first_line <- readLines(path, n = 1L)
  has_header <- grepl("barcode", first_line, ignore.case = TRUE)
  scores <- utils::read.table(path, header = has_header, stringsAsFactors = FALSE)

  label <- if (is.null(sample_id)) basename(path) else sample_id
  cells <- colnames(object)

  if (ncol(scores) >= 3L) {
    ## ---- barcode-keyed path -------------------------------------------
    colnames(scores)[1:3] <- c("barcode", "doublet_scores", "predicted_doublets")

    ## Cell names may carry a "-1" suffix, or a "-<sample>" suffix added when
    ## samples are merged. Match on the leading barcode either way.
    strip_suffix <- function(x) sub("-.*$", "", x)
    lookup <- stats::setNames(seq_len(nrow(scores)), strip_suffix(scores$barcode))
    idx <- lookup[strip_suffix(cells)]

    missing <- sum(is.na(idx))
    if (missing > 0L) {
      stop(label, ": ", missing, " of ", length(cells),
           " cells have no matching barcode in ", basename(path),
           ".\nThe score file and the count matrix are from different runs.")
    }

    object <- Seurat::AddMetaData(
      object,
      metadata = stats::setNames(scores$doublet_scores[idx], cells),
      col.name = "doublet_scores"
    )
    object <- Seurat::AddMetaData(
      object,
      metadata = stats::setNames(
        as.character(scores$predicted_doublets[idx]), cells
      ),
      col.name = "predicted_doublets"
    )
    message(label, ": joined ", length(cells), " doublet scores by barcode.")
    return(object)
  }

  ## ---- legacy positional path -----------------------------------------
  if (ncol(scores) != 2L) {
    stop(label, ": expected 2 columns (legacy) or 3 (barcoded) in ",
         basename(path), ", found ", ncol(scores), ".")
  }
  colnames(scores) <- c("doublet_scores", "predicted_doublets")

  if (nrow(scores) != length(cells)) {
    stop(label, ": doublet score file has ", nrow(scores), " rows but the ",
         "object has ", length(cells), " cells.\n",
         "This file has no barcodes, so scores can only be matched by ",
         "position, and the counts must agree exactly.\n",
         "Regenerate it with scripts/04_scrublet.py to get barcoded output.")
  }

  warning(label, ": using positional doublet join (no barcode column). ",
          "Row count matches cell count, so this is consistent, but ",
          "regenerate with scripts/04_scrublet.py to make it robust.",
          call. = FALSE)

  object <- Seurat::AddMetaData(object, metadata = scores$doublet_scores,
                                col.name = "doublet_scores")
  object <- Seurat::AddMetaData(object, metadata = scores$predicted_doublets,
                                col.name = "predicted_doublets")
  object
}

## Path to a sample's doublet score file. The pre-restructure scripts
## disagreed: the per-sample analyses read "<SAMPLE>_Doublet_scores.tsv" while
## the integrated one read a bare "Doublet_scores.tsv". Both are accepted here,
## sample-prefixed first.
doublet_score_path <- function(sample_id, cellranger_dir = CELLRANGER_DIR) {
  matrix_dir <- file.path(cellranger_dir, sample_id, "outs",
                          "filtered_feature_bc_matrix")
  candidates <- file.path(matrix_dir,
                          c(paste0(sample_id, "_Doublet_scores.tsv"),
                            "Doublet_scores.tsv"))
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0L) return(candidates[1])
  found[1]
}

## ---------------------------------------------------------------------------
## call_infected_status: label cells infected from a 2-component mixture
##
## Why this exists
## ---------------
## The pre-restructure loop assumed mixtools returned its components in a fixed
## order -- mu[1] the uninfected mode, mu[2] the high-replication mode:
##
##   ...$InfectedStatus[percent.virus > (mu[1] + sigma[1] * 2)] <- "Infected"
##   ...$InfectedStatus_groups[percent.virus > mu[2]] <- "High"
##
## normalmixEM makes no such guarantee; component order depends on
## initialization, which is random. A run that returned the components swapped
## would invert every infected/uninfected call in that sample, and nothing in
## the output would look wrong. Sorting by mu makes the labels deterministic.
##
## Also fixed here: the original fit the mixture on raw percent.virus while
## selecting which cells to fit on with a log10 cutoff, mixing the two scales.
## The fit is now done on log10 throughout, with thresholds converted back to
## the raw scale only for labelling.
## ---------------------------------------------------------------------------
call_infected_status <- function(percent_virus,
                                 log10_floor = -2.5,
                                 sd_multiplier = 2,
                                 seed = 1,
                                 label = "sample") {
  n <- length(percent_virus)
  status <- rep("Not_Infected", n)
  groups <- rep("Not_Infected", n)

  positive <- percent_virus[percent_virus > 0]
  if (length(positive) == 0L) {
    message(label, ": no viral reads in any cell; all cells Not_Infected.")
    return(list(status = status, groups = groups, fit = NULL, threshold = NA_real_))
  }

  ## Pseudocount so zero-count cells survive the log transform. The original
  ## called this "halfmin_virus" but used min(); half the minimum is the
  ## conventional choice and is what is used here.
  pseudocount <- min(positive) / 2
  log_virus <- log10(percent_virus + pseudocount)

  ## Cell selection is kept faithful to the original analysis: fit on cells
  ## whose log10 viral fraction clears log10_floor.
  ##
  ## Be aware this cutoff is pseudocount-sensitive, and was in the original
  ## too. Zero-virus cells enter the fit at log10(pseudocount), so whether they
  ## are excluded depends on how small the smallest non-zero value in the
  ## sample happens to be, not on any property of the data you meant to select
  ## on. The warning below fires when that is happening so it is at least
  ## visible. See README "Known open questions".
  fit_on <- log_virus[log_virus > log10_floor]

  if (length(fit_on) < 10L) {
    stop(label, ": only ", length(fit_on), " cells above the log10 floor of ",
         log10_floor, "; too few to fit a 2-component mixture. ",
         "Inspect this sample's percent.virus distribution by hand.")
  }

  ## A mixture needs actual spread. Without this, a fit set that is mostly one
  ## repeated value converges to two identical components with ~zero variance,
  ## and every cell above that single point gets called infected.
  if (length(unique(fit_on)) < 5L) {
    stop(label, ": fit set has only ", length(unique(fit_on)),
         " distinct values across ", length(fit_on), " cells; ",
         "a 2-component mixture is not identifiable. ",
         "Inspect this sample's percent.virus distribution by hand.")
  }

  zero_fraction <- mean(percent_virus[log_virus > log10_floor] == 0)
  if (zero_fraction > 0.5) {
    warning(sprintf(
      paste0("%s: %.0f%% of the cells entering the mixture fit have zero viral ",
             "reads -- the log10 floor of %s is not excluding them because the ",
             "pseudocount (%.3g) is large relative to it. The fitted ",
             "'uninfected' component is being driven by the pseudocount rather ",
             "than by data. Consider raising log10_floor for this sample."),
      label, 100 * zero_fraction, log10_floor, pseudocount), call. = FALSE)
  }

  set.seed(seed)  # normalmixEM initializes randomly; fix it for reproducibility
  fit <- mixtools::normalmixEM(fit_on, k = 2)

  ## Order components low -> high so mu[1]/sigma[1] is always the uninfected
  ## mode and mu[2] the high-replication mode.
  ord   <- order(fit$mu)
  mu    <- fit$mu[ord]
  sigma <- fit$sigma[ord]

  infected_cut <- mu[1] + sigma[1] * sd_multiplier
  high_cut     <- mu[2]

  status[log_virus > infected_cut] <- "Infected"

  groups[log_virus > infected_cut & log_virus <  high_cut] <- "Low"
  groups[log_virus > infected_cut & log_virus >= high_cut] <- "High"

  message(sprintf(
    "%s: mu = [%.3f, %.3f], sigma = [%.3f, %.3f] (log10 scale); infected cut %.3f; %d/%d infected",
    label, mu[1], mu[2], sigma[1], sigma[2], infected_cut,
    sum(status == "Infected"), n
  ))

  list(status = status, groups = groups, fit = fit,
       threshold = infected_cut, pseudocount = pseudocount)
}

## ---------------------------------------------------------------------------
## Apply the per-sample QC thresholds from config.R.
##
## Kept as an explicit subset() call per field rather than a string built with
## paste() so a typo is an R error, not a silently different filter.
## ---------------------------------------------------------------------------
apply_qc_filter <- function(object, sample_id,
                            thresholds = PV_QC_THRESHOLDS,
                            shared = PV_QC_SHARED) {
  t <- thresholds[[sample_id]]
  if (is.null(t)) stop("No QC thresholds defined for '", sample_id, "' in config.R")

  before <- ncol(object)
  md <- object[[]]

  keep <- rownames(md)[
    md$nCount_RNA    >  t$min_count &
    md$nCount_RNA    <  t$max_count &
    md$nFeature_RNA  >  t$min_feature &
    md$percent.ribo  <  shared$ribo_max &
    md$percent.ribo  >  t$ribo_min &
    md$percent.mt    <  shared$mt_max &
    md$percent.mt    >  shared$mt_min &
    md$doublet_scores < t$doublet_max
  ]

  object <- subset(object, cells = keep)
  message(sprintf("%s: %d -> %d cells after QC (%.1f%% kept)",
                  sample_id, before, ncol(object),
                  100 * ncol(object) / before))
  object
}

## Percentage of counts from a feature that may be absent from some samples.
## PercentageFeatureSet errors on a missing feature; this returns 0 instead.
safe_feature_percentage <- function(object, feature) {
  if (!feature %in% rownames(object)) {
    warning("Feature '", feature, "' not in the object; reporting 0%.", call. = FALSE)
    return(rep(0, ncol(object)))
  }
  result <- Seurat::PercentageFeatureSet(object, features = feature)

  ## Seurat v4 returned a one-column data.frame from PercentageFeatureSet;
  ## v5.0 returns a plain numeric vector, so the `[, 1]` this used to do fails
  ## with "incorrect number of dimensions". Accept either shape rather than
  ## pinning to one: this is the kind of return type that moves between
  ## releases, and the smoke test only caught it because it runs the real
  ## package.
  if (is.null(dim(result))) as.numeric(result) else as.numeric(result[, 1])
}
