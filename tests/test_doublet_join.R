source("analysis/helpers.R")

## Build a real Seurat object when the real package is available, so
## read_doublet_scores() is exercised against the actual AddMetaData generic.
## Under tests/run_tests.sh --stubs the stub does not export
## CreateSeuratObject, and the matrix-with-attribute mock is used instead.
##
## This matters: the matrix mock passed against the stub's plain-function
## AddMetaData but fails against the real S3 generic with "no applicable
## method for 'AddMetaData' applied to an object of class matrix" -- so a
## stub-only test could hide a real incompatibility.
USE_REAL_SEURAT <- requireNamespace("Seurat", quietly = TRUE) &&
  is.function(tryCatch(get("CreateSeuratObject", asNamespace("Seurat")),
                       error = function(e) NULL))

cat("mock objects:", if (USE_REAL_SEURAT) "real Seurat" else "matrix stub", "\n\n")

mock_object <- function(cells) {
  if (USE_REAL_SEURAT) {
    set.seed(length(cells))
    counts <- matrix(rpois(length(cells) * 6L, 20),
                     nrow = 6L,
                     dimnames = list(sprintf("GENE%d", seq_len(6L)), cells))
    return(suppressWarnings(
      Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)))
  }
  m <- matrix(0, nrow = 1, ncol = length(cells),
              dimnames = list("GENE1", cells))
  attr(m, "meta") <- list()
  m
}

md <- function(obj, field) {
  if (USE_REAL_SEURAT) return(unname(obj[[field]][, 1]))
  attr(obj, "meta")[[field]]
}
tmp <- tempfile(fileext = ".tsv")

barcodes <- sprintf("CELL%03d-1", 1:100)
scores   <- round(seq(0.01, 0.99, length.out = 100), 4)
predicted <- ifelse(scores > 0.5, "True", "False")

## --------------------------------------------------------------------------
cat("=== 1. Barcoded file, object has ALL cells ===\n")
write.table(data.frame(barcode = barcodes, doublet_score = scores,
                       predicted_doublet = predicted),
            tmp, sep = "\t", quote = FALSE, row.names = FALSE)
obj <- suppressMessages(read_doublet_scores(mock_object(barcodes), tmp, "test"))
stopifnot(identical(md(obj, "doublet_scores"), scores))
cat("PASS: all 100 scores attached in order\n\n")

## --------------------------------------------------------------------------
cat("=== 2. THE BUG: min.features dropped cell 1, scores must NOT shift ===\n")
## Seurat dropped the first barcode. Positional assignment would give every
## remaining cell the score belonging to the NEXT one.
dropped_cells <- barcodes[-1]
obj <- suppressMessages(read_doublet_scores(mock_object(dropped_cells), tmp, "test"))
got <- md(obj, "doublet_scores")
expected <- scores[-1]
shifted  <- scores[-length(scores)]   # what the old positional code produced
cat("cell CELL002 -> score", got[1], " (correct:", expected[1],
    "| old positional code would give:", shifted[1], ")\n")
stopifnot(identical(got, expected), !identical(got, shifted))
cat("PASS: barcode join is immune to the dropped cell\n\n")

## --------------------------------------------------------------------------
cat("=== 3. Cells dropped from the middle and reordered ===\n")
set.seed(7)
subset_cells <- sample(barcodes[c(-5, -20, -60, -99)])
obj <- suppressMessages(read_doublet_scores(mock_object(subset_cells), tmp, "test"))
want <- scores[match(subset_cells, barcodes)]
stopifnot(identical(md(obj, "doublet_scores"), want))
cat("PASS:", length(subset_cells), "shuffled cells all got their own score\n\n")

## --------------------------------------------------------------------------
cat("=== 4. Merged-object cell names (-<sample> suffix) ===\n")
merged_names <- sub("-1$", "-Mock_5h_PV", barcodes)
obj <- suppressMessages(read_doublet_scores(mock_object(merged_names), tmp, "test"))
stopifnot(identical(md(obj, "doublet_scores"), scores))
cat("PASS: '-Mock_5h_PV' suffix matched against '-1' in the score file\n\n")

## --------------------------------------------------------------------------
cat("=== 5. Score file from a different run is rejected ===\n")
err <- tryCatch({ read_doublet_scores(mock_object(c("OTHER-1", "CELL002-1")), tmp, "test"); "NO ERROR" },
                error = function(e) conditionMessage(e))
cat("Error:", substr(err, 1, 90), "...\n")
stopifnot(grepl("no matching barcode", err))
cat("PASS: mismatched barcodes error instead of silently dropping cells\n\n")

## --------------------------------------------------------------------------
cat("=== 6. Legacy 2-column file, counts agree -> allowed with a warning ===\n")
legacy <- tempfile(fileext = ".tsv")
write.table(data.frame(scores, predicted), legacy, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
w <- NULL
obj <- withCallingHandlers(
  suppressMessages(read_doublet_scores(mock_object(barcodes), legacy, "legacy")),
  warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
stopifnot(identical(md(obj, "doublet_scores"), scores), any(grepl("positional", w)))
cat("PASS: works, and warns that the join is positional\n\n")

## --------------------------------------------------------------------------
cat("=== 7. Legacy 2-column file, counts DISAGREE -> hard error ===\n")
cat("    (this is the silent corruption in the original code)\n")
err <- tryCatch({ suppressWarnings(read_doublet_scores(mock_object(barcodes[-1]), legacy, "legacy")); "NO ERROR" },
                error = function(e) conditionMessage(e))
cat("Error:", substr(err, 1, 95), "...\n")
stopifnot(grepl("100 rows but the object has 99 cells", err))
cat("PASS: length mismatch is fatal, not silent\n\n")

cat("=== 8. Missing file ===\n")
err <- tryCatch({ read_doublet_scores(mock_object(barcodes), "/nonexistent.tsv", "x"); "NO ERROR" },
                error = function(e) conditionMessage(e))
stopifnot(grepl("not found", err))
cat("PASS: clear error naming the scrublet step\n")
