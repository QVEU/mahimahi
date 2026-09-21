# Mirrors the real AddMetaData contract that matters here: a NAMED vector is
# matched to cells by name; an UNNAMED vector is assigned in cell order.
# The mock object is a matrix (genes x cells) so base colnames() works, with
# metadata kept in an attribute.
AddMetaData <- function(object, metadata, col.name) {
  cells <- colnames(object)
  meta <- attr(object, "meta")
  if (is.null(meta)) meta <- list()
  if (!is.null(names(metadata))) {
    meta[[col.name]] <- unname(metadata[cells])
  } else {
    if (length(metadata) != length(cells)) {
      stop("AddMetaData: length mismatch (", length(metadata), " vs ", length(cells), ")")
    }
    meta[[col.name]] <- metadata
  }
  attr(object, "meta") <- meta
  object
}

PercentageFeatureSet <- function(object, pattern = NULL, features = NULL) {
  matrix(0, nrow = ncol(object), ncol = 1)
}
