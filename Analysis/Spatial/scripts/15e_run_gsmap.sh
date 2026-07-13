#!/usr/bin/env bash
#SBATCH --job-name=gsmap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=1-15
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gsmap_run_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gsmap_run_%A_%a.err
##############################################################################
# 15e_run_gsmap.sh — Run gsMap quick_mode per sample (SLURM array).
#
# Array job: tasks 1-15 (5 GSE192741 + 10 Vu).
# Each task runs gsMap quick_mode for one spatial sample with all configured GWAS.
#
# Prerequisites:
#   - 15b: formatted .sumstats.gz + gwas_config.yaml
#   - 15c: per-sample h5ad files
#   - 15d: slice mean parquet files
#   - 15a: gsMap env + resource bundle
##############################################################################
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate gsmap

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
INPUT_DIR="${PROJECT_ROOT}/Analysis/Spatial/data/gsmap_input"
GWAS_CONFIG="${PROJECT_ROOT}/Analysis/Spatial/data/gsmap_gwas/gwas_config.yaml"
RESOURCE_DIR="${PROJECT_ROOT}/data/gsmap_resource"
WORKDIR="${PROJECT_ROOT}/Analysis/Spatial/results/gsmap"

cd "${PROJECT_ROOT}"

# ── Map SLURM_ARRAY_TASK_ID to sample ────────────────────────────────────
# Build the sample list dynamically from the h5ad files produced by 15c.
# Order: GSE192741 samples first (sorted), then Vu samples (sorted).
mapfile -t GSE_SAMPLES < <(ls "${INPUT_DIR}"/gse192741_*.h5ad 2>/dev/null | sort)
mapfile -t VU_SAMPLES  < <(ls "${INPUT_DIR}"/vu_*.h5ad 2>/dev/null | sort)

ALL_SAMPLES=("${GSE_SAMPLES[@]}" "${VU_SAMPLES[@]}")
N_TOTAL=${#ALL_SAMPLES[@]}

if [ "${N_TOTAL}" -eq 0 ]; then
    echo "ERROR: No h5ad files found in ${INPUT_DIR}. Run 15c first."
    exit 1
fi

# SLURM_ARRAY_TASK_ID is 1-based
IDX=$((SLURM_ARRAY_TASK_ID - 1))

if [ "${IDX}" -ge "${N_TOTAL}" ]; then
    echo "Task ${SLURM_ARRAY_TASK_ID} exceeds sample count (${N_TOTAL}). Exiting."
    exit 0
fi

H5AD_PATH="${ALL_SAMPLES[$IDX]}"
SAMPLE_NAME=$(basename "${H5AD_PATH}" .h5ad)

# Determine which dataset (for slice mean selection)
if [[ "${SAMPLE_NAME}" == gse192741_* ]]; then
    SLICE_MEAN="${INPUT_DIR}/gse192741_slice_mean.parquet"
    DATASET="gse192741"
else
    SLICE_MEAN="${INPUT_DIR}/vu_slice_mean.parquet"
    DATASET="vu"
fi

echo "======================================================================="
echo "  15e: gsMap Quick Mode"
echo "  Task: ${SLURM_ARRAY_TASK_ID}/${N_TOTAL}"
echo "  Sample: ${SAMPLE_NAME}"
echo "  Dataset: ${DATASET}"
echo "  H5AD: ${H5AD_PATH}"
echo "  Slice mean: ${SLICE_MEAN}"
echo "======================================================================="

# ── Validate inputs ──────────────────────────────────────────────────────
for f in "${H5AD_PATH}" "${SLICE_MEAN}" "${GWAS_CONFIG}" "${RESOURCE_DIR}/quick_mode/snp_gene_weight_matrix.h5ad"; do
    if [ ! -e "${f}" ]; then
        echo "ERROR: Required file not found: ${f}"
        exit 1
    fi
done

# ── Run gsMap quick_mode ─────────────────────────────────────────────────
# quick_mode runs all steps: find_latent_representations -> latent_to_gene
# -> generate_ldscore -> spatial_ldsc -> cauchy_combination -> report

SAMPLE_WORKDIR="${WORKDIR}/${DATASET}"
mkdir -p "${SAMPLE_WORKDIR}"

# Run gsMap; allow non-zero exit from report step (core results still valid)
gsmap quick_mode \
    --workdir "${SAMPLE_WORKDIR}" \
    --sample_name "${SAMPLE_NAME}" \
    --hdf5_path "${H5AD_PATH}" \
    --annotation condition \
    --data_layer counts \
    --gsMap_resource_dir "${RESOURCE_DIR}" \
    --sumstats_config_file "${GWAS_CONFIG}" \
    --gM_slices "${SLICE_MEAN}" \
    --max_processes 4 || {
    # Check if core results exist despite failure (report step can fail on resource limits)
    if [ -d "${SAMPLE_WORKDIR}/${SAMPLE_NAME}/cauchy_combination" ]; then
        echo "  WARNING: gsMap exited with error but core results exist — continuing"
    else
        echo "  ERROR: gsMap failed before producing core results"
        exit 1
    fi
}

echo ""
echo "======================================================================="
echo "  15e: Complete for ${SAMPLE_NAME}"
echo "======================================================================="
echo ""
echo "  Results:"
ls -lh "${SAMPLE_WORKDIR}/${SAMPLE_NAME}/" 2>/dev/null | head -20
