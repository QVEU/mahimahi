################################################################################
### SCISSORS: integrated analysis of poliovirus mutants, WT and mock
###
### Reads the per-sample Cell Ranger matrices from scripts/02_count.sh, calls
### infected cells from viral read fraction, applies per-sample QC, merges,
### regresses out cell cycle, clusters, and writes markers.
###
### Run with:  Rscript analysis/PV_mutants_integrated.R
###        or  source("analysis/PV_mutants_integrated.R")
###
### Reads are mapped against the custom reference built by scripts/00_mkref.sh
### (human GRCh38-2020-A plus the PV genome and the GFP/mRuby3 reporters).
################################################################################

library(Seurat)
library(SeuratObject)
library(clustree)
library(tidyverse)
library(ggrepel)

## Package installation is deliberately NOT done here. The pre-restructure
## script called remotes::install_version() unconditionally at the top, which
## reinstalled Seurat on every source(). Pin the environment once:
##
##   remotes::install_version("SeuratObject", "4.1.4",
##     repos = c("https://satijalab.r-universe.dev", getOption("repos")))
##   remotes::install_version("Seurat", "4.4.0",
##     repos = c("https://satijalab.r-universe.dev", getOption("repos")))
##
## This analysis targets Seurat v4. It will not run unchanged on v5, where
## layers change how merge() and the RNA assay behave.
stopifnot(packageVersion("Seurat") >= "4.0.0", packageVersion("Seurat") < "5.0.0")

## Locate this script's directory so config.R/helpers.R resolve whether the
## file is run with Rscript, sourced, or knitted.
.scissors_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
  }
  ofile <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  if (dir.exists("analysis")) return(normalizePath("analysis"))
  normalizePath(".")
}

.script_dir <- .scissors_script_dir()
source(file.path(.script_dir, "config.R"))
source(file.path(.script_dir, "helpers.R"))

################################################################################
### Load the per-sample count matrices
################################################################################

d10x.data <- sapply(PV_SAMPLE_IDS, function(id) {
  matrix_dir <- file.path(CELLRANGER_DIR, id, "outs", "filtered_feature_bc_matrix")
  if (!dir.exists(matrix_dir)) {
    stop("Count matrix not found for '", id, "': ", matrix_dir,
         "\nRun scripts/02_count.sh for this sample first.")
  }
  d10x <- Read10X(matrix_dir)
  ## Replace the 10x "-1" suffix with "-<sample>" so names.field = 2 below
  ## recovers the sample as orig.ident.
  colnames(d10x) <- paste(
    sapply(strsplit(colnames(d10x), split = "-"), "[[", 1L), id, sep = "-"
  )
  d10x
})

experiment.data <- do.call("cbind", d10x.data)

merged_seurat_object <- CreateSeuratObject(
  experiment.data,
  project      = "Mutants",
  min.cells    = 3,
  min.features = 10,
  names.field  = 2,
  names.delim  = "\\-"
)

rm(d10x.data, experiment.data)
gc(verbose = FALSE)

################################################################################
### Per-cell QC metrics, viral load and reporter expression
################################################################################

merged_seurat_object[["percent.mt"]]     <- PercentageFeatureSet(merged_seurat_object, pattern = "^MT-")
merged_seurat_object[["percent.ribo"]]   <- PercentageFeatureSet(merged_seurat_object, pattern = "^RP[SL]")
merged_seurat_object[["percent.virus"]]  <- safe_feature_percentage(merged_seurat_object, "PV")
merged_seurat_object[["mRuby"]]          <- safe_feature_percentage(merged_seurat_object, "mRuby")
merged_seurat_object[["GFP"]]            <- safe_feature_percentage(merged_seurat_object, "GFP")

################################################################################
### Call infected cells per sample from the viral read fraction
###
### Two fixes relative to the pre-restructure version:
###   - it used `break` where `next` was meant, so everything after the mock
###     sample in the list silently lost its InfectedStatus. That happened to
###     be harmless only because Mock_5h_PV is last in the vector; reordering
###     it would have broken the analysis with no error.
###   - component ordering from mixtools is now made deterministic. See
###     call_infected_status() in helpers.R.
################################################################################

sample_list <- SplitObject(merged_seurat_object, split.by = "orig.ident")

for (id in names(sample_list)) {
  if (id %in% MOCK_SAMPLE_IDS) {
    sample_list[[id]]$InfectedStatus        <- "Not_Infected"
    sample_list[[id]]$InfectedStatus_groups <- "Not_Infected"
    message(id, ": uninfected control, no mixture fit.")
    next
  }

  called <- call_infected_status(sample_list[[id]]$percent.virus, label = id)
  sample_list[[id]]$InfectedStatus        <- called$status
  sample_list[[id]]$InfectedStatus_groups <- called$groups
}

################################################################################
### Attach doublet scores and apply per-sample QC
###
### The 85 lines of copy-pasted AddMetaData blocks in the pre-restructure
### version are replaced by this loop. read_doublet_scores() joins by barcode
### where available and hard-errors on a length mismatch otherwise, instead of
### silently shifting scores onto the wrong cells. See helpers.R.
################################################################################

for (id in names(sample_list)) {
  sample_list[[id]] <- read_doublet_scores(
    sample_list[[id]],
    path      = doublet_score_path(id),
    sample_id = id
  )
  sample_list[[id]] <- apply_qc_filter(sample_list[[id]], id)
}

filtered <- merge(x = sample_list[[1]], y = sample_list[-1])

message("Merged: ", ncol(filtered), " cells across ", length(sample_list), " samples.")
print(table(filtered$predicted_doublets, filtered$orig.ident))

################################################################################
### Batch structure
###
### Mock was run alongside the CVB3 and EVA71 treatments. The IRES and MutPol
### mutants, their cotransfections and WT IRES were run together.
### Cotransfections of the transcriptional mutants were run with WT GFP.
################################################################################

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## Note: the pre-restructure scripts read
## nestorawa_forcellcycle_expressionMatrix.txt into `exp.mat` in all three
## analyses and never used it -- the cell cycle genes come from Seurat's
## built-in cc.genes. The read has been dropped.

################################################################################
### Normalize, score cell cycle, regress it out, cluster
################################################################################

x <- NormalizeData(filtered, verbose = FALSE)
x <- FindVariableFeatures(x, verbose = FALSE)
all.genes <- rownames(x)
x <- ScaleData(x, features = all.genes, verbose = FALSE)
x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes)

## Cells separate by cell cycle phase before regression.
x <- RunPCA(x, features = c(s.genes, g2m.genes), verbose = FALSE)
print(PCAPlot(x) + labs(title = "PCA on cell cycle genes, before regression"))

merged_seurat <- ScaleData(x, features = all.genes, verbose = TRUE,
                           vars.to.regress = c("S.Score", "G2M.Score"))
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(merged_seurat),
                        nfeatures.print = 10)

print(ElbowPlot(merged_seurat, ndims = 50))

merged_seurat <- FindNeighbors(merged_seurat, dims = PV_PCA_DIMS, verbose = FALSE)

################################################################################
### Choose a clustering resolution
################################################################################

clustered <- FindClusters(merged_seurat, resolution = PV_RESOLUTIONS)
print(clustree(clustered))

resolution_column <- paste0("RNA_snn_res.", PV_FINAL_RESOLUTION)
stopifnot(resolution_column %in% colnames(clustered[[]]))
Idents(clustered) <- clustered[[resolution_column]][, 1]

clustered <- RunUMAP(clustered, dims = PV_PCA_DIMS, verbose = FALSE)

print(DimPlot(clustered, group.by = resolution_column, reduction = "umap"))
print(DimPlot(clustered, group.by = resolution_column, split.by = "orig.ident",
              reduction = "umap"))

################################################################################
### Check for a batch effect between the WT GFP replicates
###
### The pre-restructure version identified the two replicate groups by
### hardcoded positional indices into Cells():
###
###   Idents(obj, cells[1561:3258])   <- 'Batch1 WT GFP'
###   Idents(obj, cells[10583:12093]) <- 'Batch2 WT GFP'
###
### Those offsets are a function of every QC threshold applied above them, so
### changing any doublet_scores cutoff would relabel the wrong cells and the
### "no significant batch effect" conclusion would be drawn from the wrong
### comparison. Selecting by orig.ident is stable under any filtering change.
###
### See the TODO on WT_REPLICATE_SAMPLE_IDS in config.R -- confirm the second
### sample before relying on this.
################################################################################

present <- intersect(WT_REPLICATE_SAMPLE_IDS, unique(clustered$orig.ident))
if (length(present) < 2L) {
  warning("Need 2 WT replicate samples for the batch check, found: ",
          paste(present, collapse = ", "), call. = FALSE)
} else {
  wt_replicates <- subset(clustered, subset = orig.ident %in% present)
  print(DimPlot(wt_replicates, group.by = "orig.ident", reduction = "umap") +
          labs(title = "WT GFP replicates, checking for batch effect"))
  print(table(wt_replicates$orig.ident, Idents(wt_replicates)))
}

print(table(clustered$orig.ident, clustered[[resolution_column]][, 1]))

################################################################################
### Save the object and metadata
################################################################################

dir.create(CELLRANGER_DIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(clustered, file = file.path(CELLRANGER_DIR, "mutsFinal.rds"))
write.csv(clustered@meta.data,
          file = file.path(CELLRANGER_DIR,
                           sprintf("scissors.metadata%s.csv", PV_FINAL_RESOLUTION)))

################################################################################
### Cluster markers
################################################################################

markers <- FindAllMarkers(clustered, only.pos = FALSE)
write.csv(markers,
          file = file.path(CELLRANGER_DIR,
                           sprintf("markers.%s.csv", PV_FINAL_RESOLUTION)))

################################################################################
### Volcano plot per cluster
################################################################################

markers$diffexpressed <- "NO"
markers$diffexpressed[markers$avg_log2FC >  1 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexpressed[markers$avg_log2FC < -1 & markers$p_val_adj < 0.05] <- "DOWN"

markers$delabel <- NA
markers$delabel[markers$diffexpressed != "NO"] <-
  markers$gene[markers$diffexpressed != "NO"]

markers$clusters <- factor(markers$cluster)

volcano <- ggplot(markers, aes(avg_log2FC, -log10(p_val_adj), shape = clusters)) +
  geom_point(size = 1, colour = "grey") +
  geom_point(data = subset(markers, diffexpressed == "UP"),   colour = "darkred") +
  geom_point(data = subset(markers, diffexpressed == "DOWN"), colour = "darkblue") +
  geom_vline(xintercept = c(-1, 1), col = "grey", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") +
  scale_shape_manual(values = seq_len(nlevels(markers$clusters))) +
  facet_wrap(~ clusters) +
  xlab(expression("log"[2] * "(Fold Change)")) +
  ylab(expression("-log"[10] * " adjusted p-value")) +
  labs(title = "Differential expression across clusters", shape = "Cluster") +
  theme_bw() +
  NoLegend()

## The pre-restructure version drew the points grey AFTER the coloured layers,
## which painted over them, and clipped the y axis at 30000 with
## coord_cartesian. Grey is now the base layer and the axis is left to the data.
print(volcano + geom_label_repel(size = 2, aes(label = delabel), max.overlaps = 1000))

################################################################################
### Marker panel
################################################################################

print(FeaturePlot(clustered, features = MARKER_PANEL))

message("Done. Object: ", file.path(CELLRANGER_DIR, "mutsFinal.rds"))
