################################################################################
### SCISSORS: replication analysis and figures
###
### Consumes results/scissors_counts.tsv.gz from the Snakemake workflow and
### produces the replication-rate figures. Replaces
### Scissors_Analysis_v4.ipynb.
###
###   Rscript analysis/scissors_replication.R
###   Rscript analysis/scissors_replication.R --counts path/to/counts.tsv.gz
###
### Computation lives in analysis/replication.R so it can be unit-tested
### (tests/test_replication_fit.R); this script is the plotting layer.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

.script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
  ofile <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  if (dir.exists("analysis")) return(normalizePath("analysis"))
  normalizePath(".")
}
SCRIPT_DIR <- .script_dir()

## Only replication.R is sourced. analysis/config.R resolves the lab share and
## stops if it cannot find it, which is right for the Seurat analyses that read
## from it -- but this script takes its input paths as arguments and never
## touches the share, so requiring it would make the workflow output
## unanalysable on a laptop.
source(file.path(SCRIPT_DIR, "replication.R"))

################################################################################
### Options
################################################################################

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}

COUNTS_FILE <- get_opt("--counts", "results/scissors_counts.tsv.gz")
FIGURE_DIR  <- get_opt("--figures", "results/figures")
## Minimum viral UMIs per cell. The notebook used UMI_count > 100 and the
## Python used CBC_readcount > 100 -- two different quantities behind the same
## number. One threshold, stated once.
MIN_UMIS    <- as.numeric(get_opt("--min-umis", "100"))
## RepliconNoInsert is the shared vector backbone: reads there cannot be
## assigned to donor or acceptor, so it is excluded from per-template fits
## (the notebook filtered it inline in each plot).
EXCLUDE     <- c("RepliconNoInsert")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
THEME <- theme_bw()

################################################################################
### Load
################################################################################

if (!file.exists(COUNTS_FILE)) {
  stop("strand counts not found: ", COUNTS_FILE,
       "\nRun the workflow first:  snakemake --cores 8")
}
counts <- read_counts_table(COUNTS_FILE)
check_counts(counts)

message(sprintf("loaded %s: %d rows, %d cells, %d samples, templates: %s",
                COUNTS_FILE, nrow(counts), length(unique(counts$CBC)),
                length(unique(counts$sample)),
                paste(sort(unique(counts$ref_name)), collapse = ", ")))

## Factor levels come from the data. The notebook hardcoded four levels in
## cell 9 ('CVB3_TT_S2_CBC.csv', 'WT_IRES_GFP_S1_CBC.csv', ...), one of which
## matched nothing it actually loaded -- factor() maps a non-match to NA, so
## that sample's point vanished from the figure instead of erroring.
sample_levels <- sort(unique(counts$sample))
counts$sample <- factor(counts$sample, levels = sample_levels)

################################################################################
### Replication estimates
###
### Two quantities, plotted on correctly-labelled axes:
###   slope     = (-)/(+)      across cells, per sample x template
###   Rep_Index = (-)/total    per cell, summarised per sample x template
################################################################################

slopes <- fit_replication_slope(counts, by = c("sample", "ref_name"),
                                min_umis = MIN_UMIS, exclude_refs = EXCLUDE)
indices <- summarise_rep_index(counts, by = c("sample", "ref_name"),
                               min_umis = MIN_UMIS, exclude_refs = EXCLUDE)

write.csv(slopes, file.path(FIGURE_DIR, "replication_slopes.csv"), row.names = FALSE)
write.csv(indices, file.path(FIGURE_DIR, "replication_indices.csv"), row.names = FALSE)

print(slopes)

fitted <- slopes[slopes$status == "ok", , drop = FALSE]
if (nrow(fitted) == 0) {
  stop("no group could be fitted. Check --min-umis (currently ", MIN_UMIS,
       ") and the status column in ", file.path(FIGURE_DIR, "replication_slopes.csv"))
}

plot_data <- counts[counts$UMI_count >= MIN_UMIS & !counts$ref_name %in% EXCLUDE, ,
                    drop = FALSE]

################################################################################
### Neg vs Pos, with the fitted slope
################################################################################

ggsave(file.path(FIGURE_DIR, "neg_vs_pos.pdf"), width = 9, height = 6,
  ggplot(plot_data, aes(Pos, Neg, colour = sample)) +
    THEME +
    geom_point(size = 0.4, alpha = 0.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    facet_wrap(~ ref_name) +
    labs(x = "Positive-strand UMIs per cell",
         y = "Negative-strand UMIs per cell",
         title = "Strand-specific UMI counts per cell",
         subtitle = sprintf("cells with >= %g viral UMIs; line slope is (-)/(+)",
                            MIN_UMIS)))

################################################################################
### Replication rate per sample
###
### The notebook labelled this axis "(-)strand/total(vRNA)" while plotting the
### Neg~Pos slope, which is (-)/(+). Both are shown here, each labelled for
### what it is. They converge only as the ratio approaches zero:
### (-)/total = r/(1+r).
################################################################################

ggsave(file.path(FIGURE_DIR, "replication_rate_slope.pdf"), width = 8, height = 4,
  ggplot(fitted, aes(x = sample, y = slope, colour = sample)) +
    THEME +
    geom_pointrange(aes(ymin = conf_low, ymax = conf_high)) +
    facet_wrap(~ ref_name) +
    ylim(0, NA) +
    labs(y = "Negative / positive strand  ((-)/(+))",
         x = NULL,
         title = "Replication rate, regression slope",
         subtitle = "point range is the 95% CI on the slope") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none"))

ggsave(file.path(FIGURE_DIR, "replication_rate_index.pdf"), width = 8, height = 4,
  ggplot(indices[indices$n_cells > 0, ], aes(x = sample, y = median_rep_index,
                                             colour = sample)) +
    THEME +
    geom_pointrange(aes(ymin = q25, ymax = q75)) +
    facet_wrap(~ ref_name) +
    ylim(0, NA) +
    labs(y = "Negative strand / total vRNA  ((-)/total)",
         x = NULL,
         title = "Replication index per cell",
         subtitle = "median with interquartile range") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none"))

################################################################################
### Distribution of the per-cell replication index
################################################################################

ggsave(file.path(FIGURE_DIR, "rep_index_distribution.pdf"), width = 8, height = 6,
  ggplot(plot_data[!is.na(plot_data$Rep_Index) & plot_data$Rep_Index > 0, ],
         aes(Rep_Index, fill = sample)) +
    THEME +
    geom_histogram(bins = 40, colour = "black", alpha = 0.5, position = "identity") +
    scale_x_log10() +
    facet_grid(sample ~ ref_name, scales = "free_y") +
    labs(x = "Negative strand / total vRNA per cell",
         y = "Cells",
         title = "Per-cell replication index") +
    theme(legend.position = "none"))

################################################################################
### Donor / acceptor comparison
###
### Requires >= 2 templates. Each cell contributes one Rep_Index per template,
### which is only meaningful because the workflow keys its wide table on
### (CBC, ref_name). The original averaged across templates, so this
### comparison was a perfect diagonal by construction.
################################################################################

templates <- setdiff(sort(unique(plot_data$ref_name)), EXCLUDE)
if (length(templates) >= 2) {
  pair <- templates[1:2]
  wide <- data.table::dcast(data.table::as.data.table(plot_data),
                            CBC + sample ~ ref_name, value.var = "Rep_Index")
  wide <- as.data.frame(wide)
  if (all(pair %in% colnames(wide))) {
    both <- wide[!is.na(wide[[pair[1]]]) & !is.na(wide[[pair[2]]]), , drop = FALSE]
    message(sprintf("%d cells carry both %s and %s", nrow(both), pair[1], pair[2]))
    if (nrow(both) > 0) {
      identical_frac <- mean(abs(both[[pair[1]]] - both[[pair[2]]]) < 1e-9)
      if (identical_frac > 0.99) {
        warning(sprintf(
          paste0("%.1f%% of co-infected cells have identical Rep_Index for %s and %s.\n",
                 "  That is the signature of the pivot_table(index=\"CBC\") bug in the ",
                 "original SCISSORS() tabulation, which averaged across templates.\n",
                 "  These counts were probably not produced by this workflow."),
          100 * identical_frac, pair[1], pair[2]), call. = FALSE)
      }
      ggsave(file.path(FIGURE_DIR, "donor_acceptor.pdf"), width = 6, height = 5.5,
        ggplot(both, aes(.data[[pair[1]]], .data[[pair[2]]], colour = sample)) +
          THEME +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
          geom_point(alpha = 0.5, size = 0.8) +
          coord_fixed() +
          labs(x = sprintf("%s replication index ((-)/total)", pair[1]),
               y = sprintf("%s replication index ((-)/total)", pair[2]),
               title = "Co-infecting template replication, per cell",
               subtitle = "dashed line is equality; points off it are template-specific"))
    }
  }
}

################################################################################
### Optional: join Seurat metadata and overlay onto the UMAP
###
### Pass --metadata=<csv> with a cell barcode column. The join reports how
### many cells matched in each direction; the notebook's inner join on a bare
### barcode string dropped non-matching cells silently, and because the
### workflow's barcodes are as-sequenced while Cell Ranger's are
### error-corrected, some loss is expected and worth seeing.
################################################################################

METADATA_FILE <- get_opt("--metadata", "")
if (nzchar(METADATA_FILE)) {
  if (!file.exists(METADATA_FILE)) stop("metadata not found: ", METADATA_FILE)
  metadata <- read_counts_table(METADATA_FILE)

  barcode_col <- get_opt("--metadata-barcode-column", colnames(metadata)[1])
  message("joining on metadata column: ", barcode_col)

  merged <- join_seurat_metadata(plot_data, metadata,
                                 counts_key = "CBC",
                                 metadata_key = barcode_col,
                                 max_unmatched_frac = 0.25)
  report <- attr(merged, "join_report")
  write.csv(data.frame(unmatched_strand_count_barcode = report$unmatched_counts),
            file.path(FIGURE_DIR, "unmatched_barcodes.csv"), row.names = FALSE)

  umap_cols <- intersect(c("UMAP_1", "UMAP_2"), colnames(merged))
  if (length(umap_cols) == 2) {
    ggsave(file.path(FIGURE_DIR, "umap_replication.pdf"), width = 10, height = 5,
      ggplot(merged, aes(UMAP_1, UMAP_2)) +
        THEME +
        geom_point(data = metadata[, c(umap_cols)], aes(UMAP_1, UMAP_2),
                   colour = "grey80", size = 0.2, inherit.aes = FALSE) +
        geom_point(aes(colour = Rep_Index, size = UMI_count), alpha = 0.8) +
        scale_colour_viridis_c(option = "magma",
                               name = "(-)/total") +
        scale_size_continuous(range = c(0.2, 2), name = "viral UMIs") +
        facet_wrap(~ ref_name) +
        coord_fixed() +
        labs(title = "Replication index on the transcriptional UMAP",
             subtitle = "grey points are all cells; coloured points carry viral UMIs"))
  } else {
    message("no UMAP_1/UMAP_2 columns in the metadata; skipping the UMAP overlay.")
  }

  cluster_col <- get_opt("--cluster-column", "seurat_clusters")
  if (cluster_col %in% colnames(merged)) {
    per_cluster <- fit_replication_slope(merged, by = c("ref_name", cluster_col),
                                         min_umis = MIN_UMIS)
    write.csv(per_cluster, file.path(FIGURE_DIR, "replication_by_cluster.csv"),
              row.names = FALSE)
    fitted_clusters <- per_cluster[per_cluster$status == "ok", , drop = FALSE]
    if (nrow(fitted_clusters) > 0) {
      ggsave(file.path(FIGURE_DIR, "replication_by_cluster.pdf"), width = 7, height = 4,
        ggplot(fitted_clusters,
               aes(x = factor(.data[[cluster_col]]), y = slope,
                   colour = factor(.data[[cluster_col]]))) +
          THEME +
          geom_pointrange(aes(ymin = conf_low, ymax = conf_high)) +
          facet_wrap(~ ref_name) +
          ylim(0, NA) +
          labs(x = "Transcriptional cluster",
               y = "Negative / positive strand  ((-)/(+))",
               title = "Replication rate by transcriptional cluster",
               subtitle = "point range is the 95% CI on the slope") +
          theme(legend.position = "none"))
    }
  }
}

message("figures and tables written to ", FIGURE_DIR)
