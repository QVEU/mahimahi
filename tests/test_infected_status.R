## Load only the mixture helper (helpers.R references Seurat:: inside other
## functions, but R does not resolve those until called).
source("analysis/helpers.R")

set.seed(42)
## Simulated sample: 2000 uninfected cells at ~0.01% viral reads, 500 infected
## cells at ~5%, 200 very high at ~30%.
percent_virus <- c(
  abs(rnorm(2000, 0.01, 0.005)),
  abs(rnorm(500,  5.0,  1.0)),
  abs(rnorm(200,  30.0, 5.0))
)
truth <- rep(c("Not_Infected", "Infected", "Infected"), c(2000, 500, 200))

cat("=== Determinism across 12 runs with a randomly-ordered mixture ===\n")
results <- replicate(12, {
  r <- suppressMessages(call_infected_status(percent_virus, label = "sim"))
  sum(r$status == "Infected")
})
cat("Infected counts per run:", paste(unique(results), collapse = ", "), "\n")
stopifnot(length(unique(results)) == 1L)
cat("PASS: identical across runs despite arbitrary component order\n\n")

cat("=== Accuracy against known truth ===\n")
r <- suppressMessages(call_infected_status(percent_virus, label = "sim"))
tab <- table(called = r$status, truth = truth)
print(tab)
acc <- sum(diag(tab[c("Infected","Not_Infected"), c("Infected","Not_Infected")])) / length(truth)
cat(sprintf("Accuracy: %.3f\n", acc))
stopifnot(acc > 0.95)
cat("PASS: infected calls match the simulated truth\n\n")

cat("=== High/Low subgrouping ===\n")
print(table(groups = r$groups, truth = truth))
stopifnot(all(c("High","Low","Not_Infected") %in% r$groups))
cat("PASS: all three groups populated\n\n")

cat("=== Mock sample with zero viral reads ===\n")
r0 <- suppressMessages(call_infected_status(rep(0, 500), label = "mock"))
stopifnot(all(r0$status == "Not_Infected"), is.null(r0$fit))
cat("PASS: all-zero input returns Not_Infected without attempting a fit\n\n")

cat("=== Too few cells above the floor ===\n")
err <- tryCatch({ suppressWarnings(call_infected_status(c(rep(0, 500), 0.9), label = "sparse")); "NO ERROR RAISED" },
                error = function(e) conditionMessage(e))
cat("Error raised:", substr(err, 1, 80), "...\n")
stopifnot(is.character(err), grepl("not identifiable|too few", err))
cat("PASS: degenerate input errors instead of returning nonsense\n")

cat("\n=== Pseudocount-swamped floor is reported, not silent ===\n")
w <- NULL
invisible(withCallingHandlers(
  suppressMessages(call_infected_status(
    c(rep(0, 500), abs(rnorm(300, 5, 1)), abs(rnorm(100, 30, 5))), label = "swamped")),
  warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") }
))
stopifnot(any(grepl("zero viral", w)))
cat("Warning raised:", substr(grep("zero viral", w, value = TRUE)[1], 1, 110), "...\n")
cat("PASS: pseudocount swamping the floor is surfaced as a warning\n")
