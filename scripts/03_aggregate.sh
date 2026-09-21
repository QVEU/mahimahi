#!/bin/bash
#SBATCH --job-name=scissors_aggr
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=logs/aggr_%j.out
#SBATCH --error=logs/aggr_%j.err
#
# Stage 3 (optional): Cell Ranger aggr across samples, with depth
# normalization. The Seurat analyses in analysis/ do NOT read this -- they read
# the per-sample filtered_feature_bc_matrix directories from stage 2 and merge
# in R. Run this only if you want Cell Ranger's own aggregated matrix.
#
# Usage:  sbatch scripts/03_aggregate.sh <aggr_csv> <run_label>
#
# Fixes relative to the pre-restructure version:
#   - "#SBATCH --ntasks-per-node = 16" had spaces around the '=', which Slurm
#     does not accept; replaced with --cpus-per-task=16.
#   - --localmem 100 was paired with #SBATCH --mem=10G, a 10x overcommit that
#     gets the job OOM-killed. Both now come from the allocation.

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

if [[ $# -lt 2 ]]; then
    echo "Usage: sbatch scripts/03_aggregate.sh <aggr_csv> <run_label>" >&2
    echo "" >&2
    echo "The aggr CSV needs a header line 'sample_id,molecule_h5' and one row" >&2
    echo "per sample pointing at that sample's outs/molecule_info.h5 from" >&2
    echo "stage 2. See README 'Aggregation CSV'." >&2
    exit 1
fi

AGGR_CSV="$1"
RUN_LABEL="$2"

[[ -f "${AGGR_CSV}" ]] || { echo "Aggregation CSV not found: ${AGGR_CSV}" >&2; exit 1; }

module load "${CELLRANGER_MODULE}"

cd "${CELLRANGER_DIR}"

cellranger aggr \
    --id="scRNAseq_aggr_${RUN_LABEL}" \
    --csv="${AGGR_CSV}" \
    --localcores="$(cellranger_localcores)" \
    --localmem="$(cellranger_localmem)"

echo "Aggregated output: ${CELLRANGER_DIR}/scRNAseq_aggr_${RUN_LABEL}/outs/"
