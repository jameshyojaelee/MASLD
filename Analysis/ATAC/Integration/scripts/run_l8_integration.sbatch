#!/bin/bash
#SBATCH --job-name=l8_integration
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/l8_integration_%j.out
#SBATCH --error=logs/l8_integration_%j.out

set -euo pipefail

# ---------------------------------------------------------------
# L8 ATAC-seq integration into the MASLD Gene Catalog (Script 35)
#
# Merges chromatin accessibility evidence from Modules 1-3 into
# the MASLD Gene Catalog as Layer 8. Handles missing modules
# gracefully — safe to run even if upstream results are partial.
#
# Usage:
#   cd Analysis/ATAC/Integration/scripts
#   sbatch run_l8_integration.sbatch              # standard run
#   sbatch run_l8_integration.sbatch --dry-run     # preview only
#   sbatch run_l8_integration.sbatch --force       # overwrite existing L8
# ---------------------------------------------------------------

echo "=== L8 ATAC-seq Integration ==="
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
SCRIPT_DIR="${PROJ_ROOT}/Analysis/ATAC/Integration/scripts"
LOG_DIR="${PROJ_ROOT}/Analysis/ATAC/Integration/logs"

mkdir -p "${LOG_DIR}"

# --- Run ---
python "${SCRIPT_DIR}/35_atac_integration.py" \
    --atlas "${PROJ_ROOT}/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv" \
    --mouse-da "${PROJ_ROOT}/Analysis/ATAC/Mouse_Bulk/results/promoter_accessibility.csv" \
    --scatac-da "${PROJ_ROOT}/Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv" \
    --chromvar "${PROJ_ROOT}/Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/chromvar_gene_annotated.csv" \
    --regulons "${PROJ_ROOT}/Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv" \
    --enhancer-links "${PROJ_ROOT}/Analysis/ATAC/Human_Multiome/scenic_plus/enhancer_gene_links.csv" \
    --ortholog-map "${PROJ_ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/human_mouse_ortholog_comparison.csv" \
    --outdir "${PROJ_ROOT}/Analysis/ATAC/Integration/results" \
    "$@"

echo ""
echo "Finished: $(date)"
