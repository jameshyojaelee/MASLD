#!/usr/bin/env bash
#SBATCH --job-name=BayesPrism
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=210G
#SBATCH --time=48:00:00
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_kallisto_%a_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_kallisto_%a_%j.err

# ---------------------------------------------------------------------------
# Stage 1: Per-cohort BayesPrism on kallisto-canonical counts
#
# Runs InstaPrism (fast BayesPrism equivalent) for each of the 5 mega-analysis
# cohorts using kallisto-derived gene-symbol count matrices.
#
# Prerequisites:
#   - Script 20_prepare_kallisto_bulk_counts.R has been run to generate
#     per-cohort counts TSVs in Analysis/Deconvolution/bulk/{COHORT}/
#   - ScaleSC reference SCE exists at reference/reference_human_scalesc_sce.rds
#
# Submit: sbatch 21_bayesprism_kallisto.sbatch
# ---------------------------------------------------------------------------

set -eo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv

set -u

export R_LIBS_USER=""
R_PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv"
export R_LIBS_SITE="${R_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${R_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
R_BIN="${R_PREFIX}/bin/R"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
REF_SCE="${PROJECT_ROOT}/Analysis/Deconvolution/reference/reference_human_scalesc_sce.rds"
BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"
RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
SCRIPTS_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/scripts"

# 5 mega-analysis cohorts (those with kallisto counts)
DATASETS=(
  "GSE126848"
  "GSE130970"
  "GSE135251"
  "GSE162694"
  "GSE213621"
)

DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"

echo "=============================="
echo "BayesPrism kallisto rerun: ${DATASET}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

COUNTS_FILE="${BULK_BASE}/${DATASET}/${DATASET}_counts.tsv"
RESULT_DIR="${RESULTS_BASE}/${DATASET}"

if [[ ! -f "${COUNTS_FILE}" ]]; then
  echo "ERROR: Missing counts file: ${COUNTS_FILE}" >&2
  echo "Run 20_prepare_kallisto_bulk_counts.R first." >&2
  exit 1
fi

if [[ ! -f "${REF_SCE}" ]]; then
  echo "ERROR: Missing SCE reference: ${REF_SCE}" >&2
  exit 1
fi

# Back up existing STAR-era BayesPrism results
PROP_OUT="${RESULT_DIR}/${DATASET}_bayesprism_proportions.tsv"
RDS_OUT="${RESULT_DIR}/${DATASET}_bayesprism_result.rds"
if [[ -f "${PROP_OUT}" ]]; then
  PROP_BACKUP="${RESULT_DIR}/${DATASET}_bayesprism_proportions_star_backup.tsv"
  if [[ ! -f "${PROP_BACKUP}" ]]; then
    cp "${PROP_OUT}" "${PROP_BACKUP}"
    echo "Backed up: ${PROP_OUT} -> $(basename ${PROP_BACKUP})"
  fi
fi
if [[ -f "${RDS_OUT}" ]]; then
  RDS_BACKUP="${RESULT_DIR}/${DATASET}_bayesprism_result_star_backup.rds"
  if [[ ! -f "${RDS_BACKUP}" ]]; then
    cp "${RDS_OUT}" "${RDS_BACKUP}"
    echo "Backed up: ${RDS_OUT} -> $(basename ${RDS_BACKUP})"
  fi
fi

echo "Running 12_run_bayesprism.R with InstaPrism..."
"${R_BIN}" --no-save --no-restore --slave \
  --file="${SCRIPTS_DIR}/12_run_bayesprism.R" \
  --args \
  "${REF_SCE}" \
  "${COUNTS_FILE}" \
  "${RESULT_DIR}" \
  "${DATASET}" \
  "celltype"

echo ""
echo "Finished: $(date)"
echo "Outputs:"
ls -lh "${RESULT_DIR}/${DATASET}_bayesprism"* 2>/dev/null || echo "  (no output files)"
