#!/bin/bash
#SBATCH --job-name=scissors_count
#SBATCH --cpus-per-task=16
#SBATCH --mem=125G
#SBATCH --time=24:00:00
#SBATCH --output=logs/count_%A_%a.out
#SBATCH --error=logs/count_%A_%a.err
#
# Stage 2: count one sample against the custom transcriptome.
#
# One sample:
#   sbatch scripts/02_count.sh Mock_5h_PV
#
# All samples in scripts/samples.tsv as an array job:
#   sbatch --array=1-$(grep -vc '^#' scripts/samples.tsv) scripts/02_count.sh
#
# The pre-restructure version hardcoded --sample, --id, --fastqs and
# --output-dir to Mock_5h and accepted a runID argument it never used, so the
# other 15 samples could not be reproduced from the repository. Sample
# parameters now come from scripts/samples.tsv.

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

# ---------------------------------------------------------------------------
# Resolve which sample to count: explicit argument, or this array task's row
# ---------------------------------------------------------------------------
sample_rows() { grep -v '^#' "${SAMPLE_TABLE}" | grep -v '^[[:space:]]*$'; }

if [[ $# -ge 1 ]]; then
    ROW="$(sample_rows | awk -F'\t' -v id="$1" '$1 == id')"
    [[ -n "${ROW}" ]] || { echo "No row for '$1' in ${SAMPLE_TABLE}" >&2; exit 1; }
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    ROW="$(sample_rows | sed -n "${SLURM_ARRAY_TASK_ID}p")"
    [[ -n "${ROW}" ]] || { echo "Array index ${SLURM_ARRAY_TASK_ID} past end of ${SAMPLE_TABLE}" >&2; exit 1; }
else
    echo "Usage: sbatch scripts/02_count.sh <output_id>" >&2
    echo "   or: sbatch --array=1-N scripts/02_count.sh" >&2
    exit 1
fi

IFS=$'\t' read -r OUTPUT_ID FASTQ_SAMPLE FASTQ_DIR <<< "${ROW}"

if [[ "${FASTQ_DIR}" == "TODO" ]]; then
    echo "${OUTPUT_ID}: fastq_dir is still TODO in ${SAMPLE_TABLE}." >&2
    echo "Fill in the stage 1 fastq_path directory for this sample first." >&2
    exit 1
fi

[[ -d "${TRANSCRIPTOME}" ]] || { echo "Transcriptome not found: ${TRANSCRIPTOME} (run 00_mkref.sh)" >&2; exit 1; }
[[ -d "${FASTQ_DIR}" ]] || { echo "FASTQ directory not found: ${FASTQ_DIR}" >&2; exit 1; }

# Confirm FASTQs for this sample prefix exist. Cell Ranger will otherwise run
# to completion on zero reads and emit an empty matrix, which is only noticed
# much later in Seurat.
if ! compgen -G "${FASTQ_DIR}/**/${FASTQ_SAMPLE}_S*_R1_*.fastq.gz" > /dev/null 2>&1 \
   && ! compgen -G "${FASTQ_DIR}/${FASTQ_SAMPLE}_S*_R1_*.fastq.gz" > /dev/null 2>&1; then
    echo "${OUTPUT_ID}: no R1 FASTQs matching prefix '${FASTQ_SAMPLE}' under ${FASTQ_DIR}" >&2
    echo "Check the fastq_sample column against SampleSheet.csv." >&2
    exit 1
fi

module load "${CELLRANGER_MODULE}"

mkdir -p "${CELLRANGER_DIR}"
cd "${CELLRANGER_DIR}"

OUTPUT_DIR="${CELLRANGER_DIR}/${OUTPUT_ID}"
if [[ -d "${OUTPUT_DIR}" ]]; then
    echo "${OUTPUT_ID}: ${OUTPUT_DIR} already exists; refusing to overwrite." >&2
    exit 1
fi

# `--create-bam` is mandatory from Cell Ranger 8.0 and rejected by 7.x, so the
# flag is derived from the installed version rather than hardcoded. An
# unparseable version fails here, before the counting starts.
if ! BAM_FLAG="$(cellranger_count_bam_flag)"; then
    exit 1
fi

echo "Sample:        ${OUTPUT_ID} (FASTQ prefix ${FASTQ_SAMPLE})"
echo "Transcriptome: ${TRANSCRIPTOME}"
echo "Cell Ranger:   $(cellranger_version)${BAM_FLAG:+  (passing ${BAM_FLAG})}"
echo "Cores/memory:  $(cellranger_localcores) cores, $(cellranger_localmem)G"

# BAM_FLAG is deliberately unquoted: it is empty on 7.x, and an empty quoted
# string would be passed as a stray argument.
# shellcheck disable=SC2086
cellranger count \
    --id="${OUTPUT_ID}" \
    --sample="${FASTQ_SAMPLE}" \
    --transcriptome="${TRANSCRIPTOME}" \
    --fastqs="${FASTQ_DIR}" \
    --localcores="$(cellranger_localcores)" \
    --localmem="$(cellranger_localmem)" \
    ${BAM_FLAG} \
    --output-dir "${OUTPUT_DIR}"

echo "Counts written to ${OUTPUT_DIR}/outs/filtered_feature_bc_matrix/"
echo "Next: scripts/04_scrublet.py ${OUTPUT_ID}"
