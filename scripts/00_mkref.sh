#!/bin/bash
#SBATCH --job-name=scissors_mkref
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/mkref_%j.out
#SBATCH --error=logs/mkref_%j.err
#
# Stage 0: build the custom transcriptome -- human GRCh38-2020-A with the
# poliovirus genome plus the GFP and mRuby3 reporters added as contigs, so
# Cell Ranger counts viral and reporter reads as ordinary features.
#
# Run once. Submit with:   sbatch scripts/00_mkref.sh
#
# The previous version of this script opened `srun --pty bash` on its first
# line, which blocks: the mkref call underneath it did not run until you
# exited that shell, and then ran on the login node rather than in the
# allocation. It only worked when the lines were pasted in by hand.

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

module load "${CELLRANGER_MODULE}"

cd "${REF_BUILD_DIR}"

for f in "${REF_FASTA}" "${REF_GTF}"; do
    [[ -f "${f}" ]] || { echo "Missing reference input: ${REF_BUILD_DIR}/${f}" >&2; exit 1; }
done

if [[ -d "${REF_GENOME_NAME}" ]]; then
    echo "Reference ${REF_BUILD_DIR}/${REF_GENOME_NAME} already exists; refusing to overwrite." >&2
    echo "Delete it first if you mean to rebuild -- downstream counts are tied to it." >&2
    exit 1
fi

cellranger mkref \
    --genome="${REF_GENOME_NAME}" \
    --fasta="${REF_FASTA}" \
    --genes="${REF_GTF}" \
    --nthreads="$(cellranger_localcores)" \
    --memgb="$(cellranger_localmem)"

echo "Reference built: ${REF_BUILD_DIR}/${REF_GENOME_NAME}"
echo "config.sh TRANSCRIPTOME resolves to: ${TRANSCRIPTOME}"
