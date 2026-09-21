#!/bin/bash
#SBATCH --job-name=scissors_mkfastq
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/mkfastq_%j.out
#SBATCH --error=logs/mkfastq_%j.err
#
# Stage 1: demultiplex the BCL run directory into per-sample FASTQs.
#
# Usage:  sbatch scripts/01_mkfastq.sh <bcl_run_dir> [run_name]
#
# Fixes relative to the pre-restructure version:
#   - "#SBATCH - cwd" removed. That is SGE's -cwd, not a Slurm option; Slurm
#     already starts in the submit directory.
#   - --localmem no longer hardcoded to 160 while #SBATCH --mem asked for 12G.
#     Both now come from the allocation (see cellranger_localmem in config.sh).
#   - --ntasks-per-node=16 replaced with --cpus-per-task=16, which is what
#     actually reserves 16 CPUs for one process.

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch scripts/01_mkfastq.sh <bcl_run_dir> [run_name]" >&2
    exit 1
fi

RUN_DIR="$1"
RUN_NAME="${2:-QVEU0056}"

[[ -d "${RUN_DIR}" ]] || { echo "BCL run directory not found: ${RUN_DIR}" >&2; exit 1; }
[[ -f "${SAMPLE_SHEET}" ]] || { echo "Sample sheet not found: ${SAMPLE_SHEET}" >&2; exit 1; }

module load "${CELLRANGER_MODULE}"
module load "${BCL2FASTQ_MODULE}"

mkdir -p "${CELLRANGER_DIR}"
cd "${CELLRANGER_DIR}"

echo "Run directory: ${RUN_DIR}"
echo "Sample sheet:  ${SAMPLE_SHEET}"
echo "Cores/memory:  $(cellranger_localcores) cores, $(cellranger_localmem)G"

cellranger mkfastq \
    --id="${RUN_NAME}" \
    --run="${RUN_DIR}" \
    --csv="${SAMPLE_SHEET}" \
    --lanes=1,2 \
    --localcores="$(cellranger_localcores)" \
    --localmem="$(cellranger_localmem)"

echo "FASTQs written under ${CELLRANGER_DIR}/${RUN_NAME}/outs/fastq_path/"
echo "Record that path in scripts/samples.tsv for stage 2."
