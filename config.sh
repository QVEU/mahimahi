#!/bin/bash
# Shared paths and Slurm helpers for the SCISSORS pipeline.
# Sourced by every script in scripts/. Edit here, not in the scripts.

# ---------------------------------------------------------------------------
# Locations on the Skyline cluster
# ---------------------------------------------------------------------------
PROJECT_ROOT="${SCISSORS_PROJECT_ROOT:-/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq}"
CELLRANGER_DIR="${PROJECT_ROOT}/CellRanger"

# Directory holding the custom reference FASTA/GTF, and the mkref output name.
REF_BUILD_DIR="${SCISSORS_REF_BUILD_DIR:-/data/lvd_qve/QVEU_Code/sequencing/template_fastas/refdata-gex-GRCh38-2020-A_PV_GFP_mRuby}"
REF_GENOME_NAME="GRCh38-2020-A_PV_GTF_mRuby"
REF_FASTA="human_pv_mrubygfp.fa"
REF_GTF="genomePVmrubygtf.gtf"

# The transcriptome 02_count.sh counts against. This is where 00_mkref.sh
# writes, so the two cannot drift apart.
#
# NOTE: the pre-restructure count script pointed somewhere else entirely --
#   .../refdata-gex-GRCh38-2020-A/fasta/GRCh38-2020-A_PV_GTF_mRuby/
# which is a different parent directory with a spurious fasta/ component, and
# is not where mkref wrote. The mkref-consistent path is used here because it
# matches the provenance comment in analysis/PV_mutants_integrated.R. If the
# published counts came from the other path, change this one line -- but the
# two references would then not be the same object, so confirm before reusing
# any existing count matrices. See README "Known open questions".
TRANSCRIPTOME="${REF_BUILD_DIR}/${REF_GENOME_NAME}"

# Illumina sample sheet consumed by 01_mkfastq.sh (see README for its format).
SAMPLE_SHEET="${CELLRANGER_DIR}/SampleSheet.csv"

# sample_id -> FASTQ directory table consumed by 02_count.sh.
SAMPLE_TABLE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts/samples.tsv"

CELLRANGER_MODULE="cellranger/7.2.0-dntehee"
BCL2FASTQ_MODULE="bcl2fastq2/2.20.0.422-orocbiu"

# ---------------------------------------------------------------------------
# Derive Cell Ranger's --localcores/--localmem from what Slurm actually granted
#
# Passing a hand-written --localmem larger than #SBATCH --mem gets the job
# OOM-killed: Cell Ranger schedules work believing it has memory the cgroup
# will not give it. Reading both back from the environment makes that class of
# mismatch impossible, so the only number to maintain is the #SBATCH line.
# ---------------------------------------------------------------------------
cellranger_localcores() {
    echo "${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-4}}"
}

cellranger_localmem() {
    # SLURM_MEM_PER_NODE is in MB. Hold back MEM_HEADROOM_GB so Cell Ranger's
    # own accounting overshoot does not cross the cgroup limit.
    local mem_mb="${SLURM_MEM_PER_NODE:-}"
    local headroom="${MEM_HEADROOM_GB:-4}"

    if [[ -z "${mem_mb}" && -n "${SLURM_MEM_PER_CPU:-}" ]]; then
        mem_mb=$(( SLURM_MEM_PER_CPU * $(cellranger_localcores) ))
    fi

    if [[ -z "${mem_mb}" ]]; then
        echo "config.sh: cannot read the Slurm memory allocation; not guessing." >&2
        echo "  Run this through sbatch, or set SLURM_MEM_PER_NODE (in MB)." >&2
        return 1
    fi

    local mem_gb=$(( mem_mb / 1024 - headroom ))
    if (( mem_gb < 1 )); then
        echo "config.sh: allocation of ${mem_mb}MB is too small for ${headroom}GB headroom." >&2
        return 1
    fi
    echo "${mem_gb}"
}

# ---------------------------------------------------------------------------
# Cell Ranger version detection
#
# The original scripts were written against Cell Ranger 7.2.0. Later releases
# changed the `count` interface: from 8.0 onwards `--create-bam` is a REQUIRED
# flag, so a 7.x command line is rejected outright by 8.x and a 8.x command
# line is rejected by 7.x. Rather than hardcode either, detect the version and
# build the flag to match, so the same script works across the range.
#
# If your cluster offers several versions, CELLRANGER_MODULE in this file
# selects one; `module avail cellranger` lists them.
# ---------------------------------------------------------------------------
cellranger_version() {
    local raw
    raw="$(cellranger --version 2>/dev/null | head -1)" || return 1
    # Formats seen across releases:
    #   "cellranger 7.2.0"
    #   "cellranger cellranger-8.0.1"
    #   "cellranger-9.0.0"
    local version
    version="$(printf '%s' "${raw}" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
    if [[ -z "${version}" ]]; then
        echo "config.sh: cannot parse a version from 'cellranger --version': ${raw}" >&2
        return 1
    fi
    printf '%s' "${version}"
}

cellranger_major() {
    local version
    version="$(cellranger_version)" || return 1
    printf '%s' "${version%%.*}"
}

# Echoes the flags `cellranger count` needs for the installed version, or
# fails loudly. Do not default to "no flag" on an unparseable version: on 8.x
# that produces a command Cell Ranger rejects after the module load, several
# minutes into a job.
cellranger_count_bam_flag() {
    local major
    if ! major="$(cellranger_major)"; then
        echo "config.sh: refusing to guess the Cell Ranger interface." >&2
        echo "  Load a cellranger module first (module load ${CELLRANGER_MODULE})." >&2
        return 1
    fi
    if (( major >= 8 )); then
        # Mandatory from 8.0. Keep the BAM: the strand-specific workflow can
        # consume a Cell Ranger BAM directly via barcode_source "tags".
        printf '%s' "--create-bam=true"
    else
        # 7.x and earlier produce a BAM unconditionally and reject the flag.
        printf '%s' ""
    fi
}
