# Extract R chunks from an .Rmd and parse each one.
args <- commandArgs(trailingOnly = TRUE)
ok <- TRUE
for (f in args) {
  lines <- readLines(f, warn = FALSE)
  starts <- grep("^```\\{r", lines)
  ends   <- grep("^```\\s*$", lines)
  n <- 0
  for (s in starts) {
    e <- ends[ends > s][1]
    if (is.na(e)) { cat(sprintf("%s: unterminated chunk at line %d\n", f, s)); ok <- FALSE; next }
    code <- lines[(s + 1):(e - 1)]
    n <- n + 1
    res <- tryCatch({ parse(text = code); NULL },
                    error = function(err) conditionMessage(err))
    if (!is.null(res)) {
      cat(sprintf("%s: chunk at line %d FAILED: %s\n", f, s, res)); ok <- FALSE
    }
  }
  cat(sprintf("%s: %d chunks, all parse OK\n", basename(f), n))
}
quit(status = if (ok) 0 else 1)
