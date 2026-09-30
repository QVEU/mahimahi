#!/usr/bin/env Rscript
## Generate a synthetic 10x share layout so the Seurat analyses can be executed.
##
##   Rscript tests/seurat/make_seurat_fixtures.R <root> [cells_per_sample]
##
## Writes, for each sample:
##   <root>/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/<id>/outs/
##       filtered_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
##       filtered_feature_bc_matrix/<id>_Doublet_scores.tsv
##
## so that SCISSORS_SHARE_ROOT=<root> makes analysis/config.R resolve to it.
##
## Written in R, not Python, so the gene universe and the count magnitudes are
## derived from the SAME definitions the analysis uses -- Seurat's cc.genes,
## and MARKER_PANEL / PV_QC_THRESHOLDS / PV_QC_SHARED from analysis/config.R.
## A Python generator would have to restate those thresholds, and could then
## drift away from the config the test is supposed to be exercising.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: make_seurat_fixtures.R <root> [cells_per_sample]")
ROOT <- normalizePath(args[[1]], mustWork = FALSE)
CELLS_PER_SAMPLE <- if (length(args) >= 2) as.integer(args[[2]]) else 80L

## RunPCA defaults to npcs = 50 and PV_PCA_DIMS asks for 23, so there has to be
## comfortably more than 50 cells per sample for the PCA to be computable.
if (CELLS_PER_SAMPLE < 60L) {
  stop("cells_per_sample must be >= 60: RunPCA(npcs = 50) needs more cells than PCs.")
}

set.seed(20260930)

dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(SCISSORS_SHARE_ROOT = ROOT)

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1]]))) else getwd()
})
REPO <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(REPO, "analysis", "config.R"))

## A subset big enough to exercise both branches of the infected-status loop:
## Mock takes the uninfected path, the other two go through the mixture fit.
PV_SAMPLES <- c("Mock_5h_PV", "WT_GFP_PV", "RFP_C109S_PV")
stopifnot(all(PV_SAMPLES %in% names(PV_QC_THRESHOLDS)))

## The two QC notebooks read a single sample each and carry their filter
## thresholds inline rather than in config.R, so those are mirrored here.
##
## NOTE: that inline placement is the same pattern config.R fixed for the PV
## samples, and these numbers now live in two places. Worth consolidating into
## config.R, but it is left alone here so this smoke test does not also change
## the notebooks' behaviour.
##   CVB3_QC.Rmd  : nCount > 7000, nFeature > 2000, percent.mt < 20, doublet < 0.53
##   EVA71_QC.Rmd : nCount > 5000, nFeature > 2000, percent.mt < 20,
##                  percent.mtribo 15-40, doublet < 0.58
NOTEBOOK_QC <- list(
  CVB3_TT  = list(min_count = 7000, max_count = Inf, min_feature = 2000,
                  ribo_min = 10, doublet_max = 0.53),
  EVA71_TT = list(min_count = 5000, max_count = Inf, min_feature = 2000,
                  ribo_min = 10, doublet_max = 0.58)
)

SAMPLES <- c(PV_SAMPLES, names(NOTEBOOK_QC))

thresholds_for <- function(sample_id) {
  if (!is.null(PV_QC_THRESHOLDS[[sample_id]])) return(PV_QC_THRESHOLDS[[sample_id]])
  if (!is.null(NOTEBOOK_QC[[sample_id]])) return(NOTEBOOK_QC[[sample_id]])
  stop("no QC thresholds known for sample '", sample_id, "'")
}

## ---------------------------------------------------------------------------
## Gene universe
##
## Every gene the analyses actually reach for has to exist or the real code
## path fails: cell-cycle scoring needs cc.genes, the RidgePlot needs four
## named genes, the closing FeaturePlot needs MARKER_PANEL, and percent.virus /
## percent.mt / percent.ribo need PV, MT- and RP[SL] features.
## ---------------------------------------------------------------------------
s_genes   <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
ridge_genes <- c("PCNA", "TOP2A", "MCM6", "MKI67")
ribo_genes <- c(sprintf("RPS%d", 2:16), sprintf("RPL%d", 2:16))
mt_genes   <- sprintf("MT-%s", c("ND1", "ND2", "CO1", "CO2", "ATP6", "CO3",
                                 "ND3", "ND4", "ND5", "ND6", "CYB", "RNR1"))
viral_genes <- c("PV", "GFP", "mRuby")

filler_n <- 6000L
filler_genes <- sprintf("GENE%04d", seq_len(filler_n))

GENES <- unique(c(s_genes, g2m_genes, ridge_genes, MARKER_PANEL,
                  ribo_genes, mt_genes, viral_genes, filler_genes))
n_genes <- length(GENES)

idx <- function(names) match(names, GENES)
RIBO_IDX  <- idx(ribo_genes)
MT_IDX    <- idx(mt_genes)
CC_IDX    <- idx(unique(c(s_genes, g2m_genes, ridge_genes)))
PANEL_IDX <- idx(MARKER_PANEL)
PV_IDX    <- idx("PV")
OTHER_IDX <- setdiff(seq_len(n_genes), c(RIBO_IDX, MT_IDX, PV_IDX))

## Three transcriptional programmes with DISJOINT marker blocks: a cell
## expresses the shared background plus its own block, and none of the other
## blocks. Without real structure every cell is statistically identical,
## Louvain returns one cluster, FindAllMarkers returns zero rows, and the smoke
## test only traverses the "no markers" guard -- never the heatmap and volcano
## it exists to cover.
##
## The separation has to be stark rather than merely present: the notebooks
## cluster at resolution 0.15-0.25, tuned for their real data, and an
## overlapping-programme fixture only separates above 0.85. Those resolutions
## are scientific parameters, so the fixture bends, not the analysis.
N_PROGRAMS <- 3L
EXCLUSIVE_PER_PROGRAM <- 800L
SHARED_N <- 2600L

programme_pool <- setdiff(OTHER_IDX, c(CC_IDX, PANEL_IDX, RIBO_IDX, MT_IDX))
stopifnot(length(programme_pool) >= SHARED_N + N_PROGRAMS * EXCLUSIVE_PER_PROGRAM)

SHARED_IDX <- programme_pool[seq_len(SHARED_N)]
EXCLUSIVE_IDX <- split(
  programme_pool[SHARED_N + seq_len(N_PROGRAMS * EXCLUSIVE_PER_PROGRAM)],
  rep(seq_len(N_PROGRAMS), each = EXCLUSIVE_PER_PROGRAM)
)

## Every cell expresses SHARED_N + EXCLUSIVE_PER_PROGRAM genes, which must
## clear the largest min_feature in play (Mock_5h_PV wants > 3000).
cat(sprintf("programmes: %d x %d exclusive genes + %d shared = %d expressed per cell\n",
            N_PROGRAMS, EXCLUSIVE_PER_PROGRAM, SHARED_N,
            SHARED_N + EXCLUSIVE_PER_PROGRAM))

cat(sprintf("gene universe: %d genes (%d cc, %d ribo, %d mt, %d panel)\n",
            n_genes, length(CC_IDX), length(RIBO_IDX), length(MT_IDX),
            length(PANEL_IDX)))

## ---------------------------------------------------------------------------
## Counts engineered to clear the REAL thresholds in config.R
##
## Read the thresholds rather than restating them, so this test fails if the
## shipped QC config and the fixtures ever stop being compatible.
## ---------------------------------------------------------------------------
make_sample_matrix <- function(sample_id, n_cells) {
  t <- thresholds_for(sample_id)
  shared <- PV_QC_SHARED

  ## Total UMIs: comfortably inside (min_count, max_count).
  target_total <- if (is.finite(t$max_count)) {
    round(mean(c(t$min_count * 1.8, min(t$max_count * 0.6, t$min_count * 3))))
  } else {
    round(t$min_count * 1.8)
  }

  ## Expressed genes are fixed by the programme structure; check rather than
  ## derive, so a threshold change in config.R fails loudly here.
  n_expressed <- length(unique(c(SHARED_IDX, EXCLUSIVE_IDX[[1]], CC_IDX, PANEL_IDX)))
  if (n_expressed <= t$min_feature) {
    stop(sprintf("%s: fixture expresses %d genes but min_feature is %g; raise SHARED_N.",
                 sample_id, n_expressed, t$min_feature))
  }

  ## Fractions inside the shared windows: ribo between t$ribo_min and 30,
  ## mt between 0.01 and 20.
  ribo_frac <- (t$ribo_min + shared$ribo_max) / 2 / 100
  mt_frac   <- 0.06

  ## Mock is the only uninfected sample; CVB3_TT and EVA71_TT are infected.
  is_mock <- sample_id %in% MOCK_SAMPLE_IDS

  i <- integer(0); j <- integer(0); x <- numeric(0)
  for (cell in seq_len(n_cells)) {
    total <- round(target_total * runif(1, 0.9, 1.1))

    n_ribo <- round(total * ribo_frac)
    n_mt   <- round(total * mt_frac)

    ## Viral load: mock has none. Infected samples get a bimodal
    ## distribution so mixtools has two components to find, with about a
    ## quarter of cells uninfected -- below the 50% that would trip the
    ## pseudocount warning in call_infected_status().
    n_virus <- if (is_mock) 0L else {
      r <- runif(1)
      if (r < 0.25) 0L
      else if (r < 0.65) round(total * runif(1, 0.002, 0.01))
      else round(total * runif(1, 0.04, 0.12))
    }

    n_rest <- total - n_ribo - n_mt - n_virus
    stopifnot(n_rest > n_expressed)

    ## This cell's programme, cycling deterministically so every sample holds
    ## all three in equal proportion.
    program <- ((cell - 1L) %% N_PROGRAMS) + 1L
    exclusive <- EXCLUSIVE_IDX[[program]]

    ## Shared background plus this programme's exclusive block, and nothing
    ## from the other two blocks -- that zero is what makes the profiles
    ## separable at low resolution.
    expressed <- c(SHARED_IDX, exclusive, CC_IDX, PANEL_IDX)
    expressed <- unique(expressed)
    counts <- rep(1L, length(expressed))
    extra <- n_rest - length(expressed)
    stopifnot(extra > 0)

    ## Load the remainder onto the exclusive block, so those genes are both
    ## present-only-here and highly expressed, plus the cell-cycle and marker
    ## genes so CellCycleScoring and the closing FeaturePlot have real signal.
    weighted <- c(rep(exclusive, 12L),
                  rep(idx(unique(c(s_genes, g2m_genes))), 2L),
                  PANEL_IDX,
                  sample(SHARED_IDX, min(400L, length(SHARED_IDX))))
    hits <- table(sample(weighted, extra, replace = TRUE))
    pos <- match(as.integer(names(hits)), expressed)
    counts[pos] <- counts[pos] + as.integer(hits)

    gi <- expressed; gx <- counts
    if (n_ribo > 0) {
      h <- table(sample(RIBO_IDX, n_ribo, replace = TRUE))
      gi <- c(gi, as.integer(names(h))); gx <- c(gx, as.integer(h))
    }
    if (n_mt > 0) {
      h <- table(sample(MT_IDX, n_mt, replace = TRUE))
      gi <- c(gi, as.integer(names(h))); gx <- c(gx, as.integer(h))
    }
    if (n_virus > 0) {
      gi <- c(gi, PV_IDX); gx <- c(gx, n_virus)
    }

    i <- c(i, gi); j <- c(j, rep(cell, length(gi))); x <- c(x, gx)
  }

  m <- sparseMatrix(i = i, j = j, x = x, dims = c(n_genes, n_cells))
  rownames(m) <- GENES
  colnames(m) <- sprintf("%s-1", vapply(seq_len(n_cells), function(k)
    paste0(sample(c("A", "C", "G", "T"), 16, replace = TRUE), collapse = ""),
    character(1)))
  m
}

write_10x <- function(m, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  mtx <- file.path(dir, "matrix.mtx")
  Matrix::writeMM(m, mtx)
  system2("gzip", c("-f", shQuote(mtx)))

  ## Cell Ranger v3 features.tsv.gz: id, symbol, type. Read10X reads column 2.
  feat <- gzfile(file.path(dir, "features.tsv.gz"), "wt")
  writeLines(paste(sprintf("ENSGSYN%08d", seq_len(nrow(m))), rownames(m),
                   "Gene Expression", sep = "\t"), feat)
  close(feat)

  bc <- gzfile(file.path(dir, "barcodes.tsv.gz"), "wt")
  writeLines(colnames(m), bc)
  close(bc)
}

write_doublet_scores <- function(m, dir, sample_id) {
  t <- thresholds_for(sample_id)
  ## Comfortably below this sample's cutoff so QC keeps every cell; the
  ## doublet filter itself is unit-tested in tests/test_doublet_join.R.
  scores <- runif(ncol(m), 0.05, t$doublet_max * 0.6)
  out <- data.frame(barcode = colnames(m),
                    doublet_score = round(scores, 6),
                    predicted_doublet = ifelse(scores > 0.5, "True", "False"),
                    stringsAsFactors = FALSE)
  write.table(out, file.path(dir, sprintf("%s_Doublet_scores.tsv", sample_id)),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

for (sample_id in SAMPLES) {
  m <- make_sample_matrix(sample_id, CELLS_PER_SAMPLE)
  dir <- file.path(CELLRANGER_DIR, sample_id, "outs", "filtered_feature_bc_matrix")
  write_10x(m, dir)
  write_doublet_scores(m, dir, sample_id)

  ## Report the QC metrics against the thresholds this sample will be filtered
  ## on, so a fixture that cannot survive QC is obvious here rather than as an
  ## empty object three steps later.
  t <- thresholds_for(sample_id)
  total <- Matrix::colSums(m)
  feats <- Matrix::colSums(m > 0)
  ribo <- 100 * Matrix::colSums(m[RIBO_IDX, , drop = FALSE]) / total
  mt   <- 100 * Matrix::colSums(m[MT_IDX, , drop = FALSE]) / total
  virus <- 100 * m[PV_IDX, ] / total
  cat(sprintf(
    paste0("%-14s %3d cells  nCount %6.0f-%6.0f (>%g)  nFeature %4.0f-%4.0f (>%g)  ",
           "ribo %.1f-%.1f (%g-%g)  mt %.1f-%.1f  virus %.2f-%.2f%%\n"),
    sample_id, ncol(m), min(total), max(total), t$min_count,
    min(feats), max(feats), t$min_feature,
    min(ribo), max(ribo), t$ribo_min, PV_QC_SHARED$ribo_max,
    min(mt), max(mt), min(virus), max(virus)))

  stopifnot(all(total > t$min_count), all(total < t$max_count),
            all(feats > t$min_feature),
            all(ribo > t$ribo_min), all(ribo < PV_QC_SHARED$ribo_max),
            all(mt > PV_QC_SHARED$mt_min), all(mt < PV_QC_SHARED$mt_max))
}

cat(sprintf("\nfixtures written under %s\n", CELLRANGER_DIR))
cat(sprintf("run with: SCISSORS_SHARE_ROOT=%s Rscript analysis/PV_mutants_integrated.R\n", ROOT))
