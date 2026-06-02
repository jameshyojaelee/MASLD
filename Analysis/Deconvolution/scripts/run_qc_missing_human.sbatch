#!/usr/bin/env bash
#SBATCH --job-name=qc_missing
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/qc_missing_%a_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/qc_missing_%a_%j.err

# QC plots for the 3 human datasets after BayesPrism backfill.
# Runs 08_qc_music.R (MuSiC proportions + metadata -> heatmap, bar plots, etc.)
# Resources: 32G / 2h provides ample margin for ggplot2 plotting workloads.
#
# DEPENDENCY: Runs after run_bayesprism_missing_human.sbatch (dependency chain).

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV_PREFIX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv

BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"
RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
SCRIPTS_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/scripts"

R_BIN="${ENV_PREFIX}/bin/R"
export R_LIBS_USER=""
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

DATASETS=(
  "GSE126848"
  "PRJNA512027"
  "GSE167523"
)

DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"

echo "=============================="
echo "QC plots: ${DATASET}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

MUSIC_PROP="${RESULTS_BASE}/${DATASET}/${DATASET}_music_prop_weighted.tsv"
METADATA="${BULK_BASE}/${DATASET}/${DATASET}_metadata.tsv"
QC_DIR="${RESULTS_BASE}/${DATASET}/qc"

if [[ ! -f "${MUSIC_PROP}" ]]; then
  echo "ERROR: Missing MuSiC proportions: ${MUSIC_PROP}" >&2
  exit 1
fi

if [[ ! -f "${METADATA}" ]]; then
  echo "ERROR: Missing metadata: ${METADATA}" >&2
  exit 1
fi

mkdir -p "${QC_DIR}"

echo "Running 08_qc_music.R..."
"${R_BIN}" --no-save --no-restore --slave \
  --file="${SCRIPTS_DIR}/08_qc_music.R" \
  --args \
  "${MUSIC_PROP}" \
  "${METADATA}" \
  "${QC_DIR}" \
  "${DATASET}"

echo "Finished: $(date)"
echo "QC outputs in ${QC_DIR}:"
ls -lh "${QC_DIR}/" 2>/dev/null || echo "  (no output files found)"
