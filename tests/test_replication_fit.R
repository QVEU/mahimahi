## Tests for analysis/replication.R -- base R only, no external packages.
source("analysis/replication.R")

fails <- 0
ok <- function(label, cond) {
  cat(if (isTRUE(cond)) "  PASS  " else "  FAIL  ", label, "\n", sep = "")
  if (!isTRUE(cond)) fails <<- fails + 1
}

## --------------------------------------------------------------------------
cat("=== 0. The original fitSet() shape: does `break` in a function work? ===\n")
original_shape <- function(x) {
  first <- 1
  break                      # exactly where fitSet() had it
  second <- 2
  return(second)
}
err <- tryCatch({ original_shape(1); "NO ERROR" },
                error = function(e) conditionMessage(e))
cat("  R says:", err, "\n")
ok("bare `break` in a function is an error (so fitSet() never returned)",
   grepl("no loop for break", err))

## --------------------------------------------------------------------------
cat("\n=== 1. Recover a known slope ===\n")
set.seed(11)
make_group <- function(sample, ref, n, true_slope, noise = 2) {
  pos <- round(runif(n, 50, 2000))
  data.frame(CBC = sprintf("CELL%04d", seq_len(n)), ref_name = ref, sample = sample,
             Pos = pos,
             Neg = pmax(0, round(true_slope * pos + rnorm(n, 0, noise))),
             UMI_count = pos, stringsAsFactors = FALSE)
}
counts <- rbind(make_group("S1", "eGFP", 200, 0.05),
                make_group("S1", "mRuby3", 200, 0.20))
counts$Rep_Index <- counts$Neg / (counts$Pos + counts$Neg)

fits <- fit_replication_slope(counts)
print(fits[, c("sample", "ref_name", "n_cells", "slope", "conf_low", "conf_high", "status")])
egfp <- fits[fits$ref_name == "eGFP", ]
ruby <- fits[fits$ref_name == "mRuby3", ]
ok("returns a data frame with a slope (the thing fitSet could not do)",
   is.data.frame(fits) && all(c("slope", "slope_se") %in% colnames(fits)))
ok("eGFP slope within 0.005 of the true 0.05", abs(egfp$slope - 0.05) < 0.005)
ok("mRuby3 slope within 0.005 of the true 0.20", abs(ruby$slope - 0.20) < 0.005)
ok("true slope inside the 95% CI (eGFP)", egfp$conf_low < 0.05 && egfp$conf_high > 0.05)
ok("true slope inside the 95% CI (mRuby3)", ruby$conf_low < 0.20 && ruby$conf_high > 0.20)

## --------------------------------------------------------------------------
cat("\n=== 2. CI is wider than +/-1 SE ===\n")
cat(sprintf("  slope %.5f  se %.5f  +/-1SE = [%.5f, %.5f]  95%%CI = [%.5f, %.5f]\n",
            egfp$slope, egfp$slope_se, egfp$slope - egfp$slope_se,
            egfp$slope + egfp$slope_se, egfp$conf_low, egfp$conf_high))
ok("95% CI is wider than the +/-1 SE the notebook plotted",
   (egfp$conf_high - egfp$conf_low) > 2 * egfp$slope_se)

## --------------------------------------------------------------------------
cat("\n=== 3. Degenerate groups are reported, not fatal ===\n")
tiny <- rbind(counts,
              data.frame(CBC = c("X1", "X2"), ref_name = "TinyRef", sample = "S1",
                         Pos = c(10, 20), Neg = c(1, 2), UMI_count = c(10, 20),
                         Rep_Index = c(0.09, 0.09), stringsAsFactors = FALSE),
              data.frame(CBC = paste0("Z", 1:8), ref_name = "FlatRef", sample = "S1",
                         Pos = rep(100, 8), Neg = 0:7, UMI_count = rep(100, 8),
                         Rep_Index = 0, stringsAsFactors = FALSE),
              data.frame(CBC = paste0("W", 1:8), ref_name = "NoNegRef", sample = "S1",
                         Pos = seq(10, 80, 10), Neg = rep(0, 8), UMI_count = seq(10, 80, 10),
                         Rep_Index = 0, stringsAsFactors = FALSE))
res <- suppressMessages(fit_replication_slope(tiny))
print(res[, c("ref_name", "n_cells", "slope", "status")])
ok("small group flagged insufficient_cells, not an error",
   res$status[res$ref_name == "TinyRef"] == "insufficient_cells")
ok("zero-variance Pos flagged no_variance_in_Pos",
   res$status[res$ref_name == "FlatRef"] == "no_variance_in_Pos")
ok("all-zero Neg flagged no_negative_strand",
   res$status[res$ref_name == "NoNegRef"] == "no_negative_strand")
ok("the fittable groups still succeeded",
   all(res$status[res$ref_name %in% c("eGFP", "mRuby3")] == "ok"))

## --------------------------------------------------------------------------
cat("\n=== 4. Coefficients extracted by name, not position ===\n")
## A 1-row-per-group fit changes the coefficient matrix shape. Positional
## indexing ([2], [4]) would read the wrong cell; named indexing cannot.
coefs <- coef(summary(lm(Neg ~ Pos, data = counts[counts$ref_name == "eGFP", ])))
ok("positional [2] happens to equal the named slope here",
   abs(coefs[2] - coefs["Pos", "Estimate"]) < 1e-12)
ok("positional [4] happens to equal the named slope SE here",
   abs(coefs[4] - coefs["Pos", "Std. Error"]) < 1e-12)
no_intercept <- coef(summary(lm(Neg ~ Pos + 0, data = counts[counts$ref_name == "eGFP", ])))
cat("  with no intercept the matrix is", nrow(no_intercept), "x", ncol(no_intercept),
    "and [2] is now", sprintf("%.5f", no_intercept[2]),
    "while the slope is", sprintf("%.5f", no_intercept["Pos", "Estimate"]), "\n")
ok("positional [2] is WRONG once the matrix shape changes",
   abs(no_intercept[2] - no_intercept["Pos", "Estimate"]) > 1e-9)

## --------------------------------------------------------------------------
cat("\n=== 5. slope is (-)/(+), Rep_Index is (-)/total -- different numbers ===\n")
idx <- summarise_rep_index(counts)
print(idx)
for (ref in c("eGFP", "mRuby3")) {
  r <- fits$slope[fits$ref_name == ref]
  expected_total <- r / (1 + r)
  observed <- idx$mean_rep_index[idx$ref_name == ref]
  cat(sprintf("  %-7s slope (-)/(+) = %.4f ; r/(1+r) = %.4f ; mean Rep_Index = %.4f\n",
              ref, r, expected_total, observed))
}
r_ruby <- fits$slope[fits$ref_name == "mRuby3"]
ok("at r=0.2 the two quantities differ by >10% (so the axis label matters)",
   abs(r_ruby - r_ruby / (1 + r_ruby)) / r_ruby > 0.10)

## --------------------------------------------------------------------------
cat("\n=== 6. Missing columns fail with a clear message ===\n")
err <- tryCatch({ fit_replication_slope(data.frame(a = 1)); "NO ERROR" },
                error = function(e) conditionMessage(e))
cat("  ", substr(err, 1, 78), "...\n", sep = "")
ok("missing columns named explicitly", grepl("missing required column", err))

## --------------------------------------------------------------------------
cat("\n=== 7. Seurat join reports losses in both directions ===\n")
sc_counts <- data.frame(CBC = c("AAAA", "BBBB", "CCCC", "DDDD"),
                        ref_name = "eGFP", sample = "S1",
                        Pos = c(100, 200, 300, 400), Neg = c(5, 10, 15, 20),
                        Rep_Index = 0.05, UMI_count = c(105, 210, 315, 420),
                        stringsAsFactors = FALSE)
meta <- data.frame(cell_barcode = c("AAAA-1", "BBBB-1", "EEEE-1"),
                   seurat_clusters = c(0, 1, 2), stringsAsFactors = FALSE)
msgs <- capture.output(
  joined <- join_seurat_metadata(sc_counts, meta, max_unmatched_frac = 0.9),
  type = "message")
cat("  ", paste(msgs, collapse = "\n  "), "\n", sep = "")
rep <- attr(joined, "join_report")
ok("matched the 2 shared cells", rep$matched_cells == 2)
ok("'-1' suffix handled", all(c("AAAA", "BBBB") %in% strip_barcode_suffix(joined$CBC)))
ok("reports the 2 unmatched strand-count cells", length(rep$unmatched_counts) == 2)
ok("reports the 1 Seurat cell with no counts", length(rep$unmatched_metadata) == 1)

cat("\n=== 8. Excess unmatched fraction can be made fatal ===\n")
err <- tryCatch({
  suppressMessages(join_seurat_metadata(sc_counts, meta, max_unmatched_frac = 0.1,
                                        on_excess = "error"))
  "NO ERROR"
}, error = function(e) conditionMessage(e))
ok("errors when loss exceeds the threshold", grepl("did not match", err))
ok("names example unmatched barcodes", grepl("CCCC|DDDD", err))

cat("\n")
if (fails > 0) { cat("FAILED:", fails, "\n"); quit(status = 1) }
cat("All replication fit tests passed.\n")
