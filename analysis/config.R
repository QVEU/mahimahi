## Shared paths for the SCISSORS R analyses.
## Sourced by every script in analysis/. Edit here, not in the scripts.

## ---------------------------------------------------------------------------
## Resolve the lab share root
##
## The pre-restructure scripts referred to the same share four different ways:
##   /Volumes/lvd_qve   (10x)   macOS SMB mount
##   /Volumes/LVD_QVE   (2x)    same mount, different case -- broke on any
##                              case-sensitive volume
##   /data/lvd_qve      (4x)    Skyline cluster
##   ~/lab_share        (2x)    another local mount
## One script read from ~/lab_share and wrote to /data/lvd_qve in the same run.
## Resolving once here means an analysis moves between the cluster and a laptop
## without editing any paths.
## ---------------------------------------------------------------------------
resolve_share_root <- function() {
  override <- Sys.getenv("SCISSORS_SHARE_ROOT", unset = NA)
  if (!is.na(override) && nzchar(override)) {
    if (!dir.exists(override)) {
      stop("SCISSORS_SHARE_ROOT is set to a path that does not exist: ", override)
    }
    return(normalizePath(override))
  }

  candidates <- c(
    "/data/lvd_qve",          # Skyline cluster
    "~/lab_share",            # local mount
    "/Volumes/lvd_qve",       # macOS SMB mount
    "/Volumes/LVD_QVE"        # same mount on a case-insensitive volume
  )
  for (path in candidates) {
    expanded <- path.expand(path)
    if (dir.exists(expanded)) return(normalizePath(expanded))
  }

  stop("Could not find the lab share. Tried:\n  ",
       paste(candidates, collapse = "\n  "),
       "\nMount it, or set SCISSORS_SHARE_ROOT to its location.")
}

SHARE_ROOT     <- resolve_share_root()
PROJECT_DIR    <- file.path(SHARE_ROOT, "Projects", "PTD_StrandSpecificCounting_scRNAseq")
CELLRANGER_DIR <- file.path(PROJECT_DIR, "CellRanger")
RESULTS_DIR    <- file.path(SHARE_ROOT, "Projects", "CM_kb")

## ---------------------------------------------------------------------------
## Sample sets
## ---------------------------------------------------------------------------

## Poliovirus samples for the integrated analysis. Must match the output_id
## column of scripts/samples.tsv. Mock is kept last only for readability now --
## PV_mutants_integrated.R no longer depends on its position (see NEWS).
PV_SAMPLE_IDS <- c(
  "RFP_C109S_PV",
  "WT_GFP_PV",
  "WT_GFP_RFP_Y88P_PV",
  "WT_GFP_RFP_D177A_PV",
  "RFP_Y88P_PV",
  "RFP_D177A_PV",
  "WT_GFP_RFP_C109S_PV",
  "WT_IRES_GFP_PV",
  "Del_IRES_mRuby3_PV",
  "WT_IRES_GFP_Del_IRES_mRuby3_PV",
  "WT_IRES_mRuby3_MutPol_PV",
  "Del_IRES_mRuby3_MutPol_PV",
  "WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV",
  "Mock_5h_PV"
)

## Uninfected controls, excluded from viral-load mixture modelling.
MOCK_SAMPLE_IDS <- c("Mock_5h_PV")

## Samples carrying a WT GFP virus, compared against each other to check for a
## batch effect.
##
## TODO: confirm this list. The pre-restructure script identified the two
## groups by hardcoded cell indices -- merged cells [1561:3258] as "Batch1 WT
## GFP" and [10583:12093] as "Batch2 WT GFP". Given the merge order,
## [1561:3258] (1698 cells) is WT_GFP_PV. The second range (1511 cells) falls
## around the eighth sample, most plausibly WT_IRES_GFP_PV, but that cannot be
## confirmed without the per-sample cell counts from that run. Set this to the
## two samples you actually intended to compare.
WT_REPLICATE_SAMPLE_IDS <- c("WT_GFP_PV", "WT_IRES_GFP_PV")

## Per-sample QC thresholds, carried over verbatim from the pre-restructure
## scripts. These were hand-tuned against each sample's violin plots and
## doublet score distribution, so they are data, not defaults -- do not
## "tidy" them into a single shared threshold.
##
## doublet_max applies to the Scrublet score. If you regenerate scores with
## scripts/04_scrublet.py, re-tune these against the new distributions.
PV_QC_THRESHOLDS <- list(
  Mock_5h_PV                           = list(min_count = 5000, max_count = 40000, min_feature = 3000, ribo_min = 10, doublet_max = 0.57),
  RFP_C109S_PV                         = list(min_count = 5000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.51),
  WT_GFP_PV                            = list(min_count = 5000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.54),
  WT_GFP_RFP_Y88P_PV                   = list(min_count = 4000, max_count = Inf,   min_feature = 1500, ribo_min = 5,  doublet_max = 0.55),
  WT_GFP_RFP_D177A_PV                  = list(min_count = 4000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.50),
  RFP_Y88P_PV                          = list(min_count = 3000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.47),
  RFP_D177A_PV                         = list(min_count = 2000, max_count = Inf,   min_feature = 1000, ribo_min = 5,  doublet_max = 0.52),
  WT_GFP_RFP_C109S_PV                  = list(min_count = 2000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.46),
  WT_IRES_GFP_PV                       = list(min_count = 3000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.45),
  Del_IRES_mRuby3_PV                   = list(min_count = 3000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.45),
  WT_IRES_GFP_Del_IRES_mRuby3_PV       = list(min_count = 3000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.53),
  WT_IRES_mRuby3_MutPol_PV             = list(min_count = 2000, max_count = Inf,   min_feature = 1000, ribo_min = 5,  doublet_max = 0.57),
  Del_IRES_mRuby3_MutPol_PV            = list(min_count = 2000, max_count = Inf,   min_feature = 1000, ribo_min = 5,  doublet_max = 0.57),
  WT_IRES_mRuby3_MutPol_WT_IRES_GFP_PV = list(min_count = 5000, max_count = Inf,   min_feature = 2000, ribo_min = 5,  doublet_max = 0.54)
)

## Thresholds shared by every PV sample.
PV_QC_SHARED <- list(ribo_max = 30, mt_max = 20, mt_min = 0.01)

## Clustering parameters for the integrated analysis.
PV_PCA_DIMS       <- 1:23
PV_RESOLUTIONS    <- seq(0.1, 1, by = 0.1)
PV_FINAL_RESOLUTION <- 0.5

## Interferon-stimulated and stress-response genes plotted at the end of each
## analysis.
MARKER_PANEL <- c("ISG15", "OASL", "CCN1", "WFDC1", "AKR1C2", "IFIT3",
                  "CRLF1", "AKR1C1", "HIST1H2AE", "IFIT2", "UBE2C",
                  "CD74", "MT1G", "HIST1H4C")
