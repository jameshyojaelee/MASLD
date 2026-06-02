#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/scenic_grn_%j.out
#SBATCH --error=logs/scenic_grn_%j.out

set -euo pipefail

# ---------------------------------------------------------------
# Module 3: Gene-Activity-Based Hepatocyte GRN
#
# Builds TF regulons and enhancer-gene links using SnapATAC2 gene
# activity matrix as RNA proxy (separate ATAC/RNA assays cannot be
# barcode-matched, so native SCENIC+ is not feasible).
#
# Outputs 4 CSVs to scenic_plus/ for L8 integration (Script 35).
#
# Usage:
#   cd Analysis/ATAC/Human_Multiome/scripts
#   sbatch run_scenic_grn.sbatch
# ---------------------------------------------------------------

echo "=== Module 3: Gene-Activity-Based Hepatocyte GRN ==="
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Node:      $(hostname)"
echo "CPUs:      ${SLURM_CPUS_PER_TASK}"
echo "Memory:    ${SLURM_MEM_PER_NODE}MB"
echo "Started:   $(date)"
echo ""

# --- Environment ---
module purge 2>/dev/null || true
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# --- Paths ---
PROJ_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
MULTIOME_DIR="${PROJ_ROOT}/Analysis/ATAC/Human_Multiome"
SCRIPT_DIR="${MULTIOME_DIR}/scripts"
RESULTS_DIR="${MULTIOME_DIR}/results/snapatac2"
OUTPUT_DIR="${MULTIOME_DIR}/scenic_plus"
LOG_DIR="${MULTIOME_DIR}/logs"

GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"

mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}"

# Prefer fixed h5ad (unique obs_names) if available
ATAC_H5AD="${RESULTS_DIR}/snapatac2_processed_fixed.h5ad"
if [[ ! -f "${ATAC_H5AD}" ]]; then
    ATAC_H5AD="${RESULTS_DIR}/snapatac2_processed.h5ad"
fi
echo "Using ATAC h5ad: ${ATAC_H5AD}"

# --- Run ---
python "${SCRIPT_DIR}/04b_scenic_grn_from_activity.py" \
    --atac-h5ad "${ATAC_H5AD}" \
    --anndataset "${RESULTS_DIR}/combined.h5ads" \
    --gtf "${GTF}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-cells 15000 \
    --window 500000 \
    --min-corr 0.1 \
    --top-peaks 50000

echo ""
echo "=== Verifying outputs ==="
for f in hepatocyte_regulons.csv enhancer_gene_links.csv regulon_activity_scores.csv disease_regulons.csv; do
    if [[ -f "${OUTPUT_DIR}/${f}" ]]; then
        lines=$(wc -l < "${OUTPUT_DIR}/${f}")
        echo "  OK: ${f} (${lines} lines)"
    else
        echo "  MISSING: ${f}"
    fi
done

echo ""
echo "Finished: $(date)"
