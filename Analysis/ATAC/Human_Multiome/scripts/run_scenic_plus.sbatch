#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=logs/scenic_plus_%j.out
#SBATCH --error=logs/scenic_plus_%j.err

# =============================================================================
# SCENIC+ Hepatocyte Gene Regulatory Network (Module 3)
#
# Runs the full SCENIC+ pipeline on paired snRNA + snATAC multiome data
# (GSE244832) to build a hepatocyte-specific GRN.
#
# Also runs chromVAR TF motif analysis (Module 2c) as a prerequisite.
#
# Prerequisites:
#   - SnapATAC2 processed AnnData: results/snapatac2/snapatac2_processed.h5ad
#   - CellRanger ARC outputs in cellranger_arc/{donor}/outs/
#   - atac_env conda environment with scenicplus, pycisTopic, pycistarget
#
# Usage:
#   cd Analysis/ATAC/Human_Multiome
#   mkdir -p logs
#   sbatch scripts/run_scenic_plus.sbatch
#
# Estimated runtime: 12-24h on GPU with 15K cells
# =============================================================================

set -euo pipefail

# --- Configuration ----------------------------------------------------------
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WORK_DIR="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome"
SCRIPT_DIR="${WORK_DIR}/scripts"

ATAC_INPUT="${WORK_DIR}/results/snapatac2/snapatac2_processed.h5ad"
RNA_DIR="${WORK_DIR}/cellranger_arc"
OUTPUT_DIR="${WORK_DIR}/results/scenic_plus"
CHROMVAR_DIR="${WORK_DIR}/results/chromvar"

N_CELLS=15000
N_TOPICS="10,20,30,40"
N_ITER=500
N_HVG=3000

# --- Environment setup ------------------------------------------------------
echo "============================================================"
echo "SCENIC+ Hepatocyte GRN Pipeline"
echo "  Job ID:     ${SLURM_JOB_ID}"
echo "  Node:       $(hostname)"
echo "  GPUs:       ${CUDA_VISIBLE_DEVICES:-none}"
echo "  CPUs:       ${SLURM_CPUS_PER_TASK}"
echo "  Memory:     ${SLURM_MEM_PER_NODE:-256G}"
echo "  Started:    $(date)"
echo "============================================================"

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate atac_env

# Fix LD_LIBRARY_PATH for pip-installed NVIDIA libs (if present)
CONDA_PREFIX_PATH="$(micromamba info | grep 'env location' | awk '{print $NF}')"
if [[ -d "${CONDA_PREFIX_PATH}/lib" ]]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX_PATH}/lib:${LD_LIBRARY_PATH:-}"
fi

# Report GPU status
echo ""
echo "GPU status:"
nvidia-smi || echo "  nvidia-smi not available"
echo ""

python -c "
import torch
if torch.cuda.is_available():
    print(f'  PyTorch CUDA: {torch.cuda.get_device_name(0)}')
    print(f'  CUDA version: {torch.version.cuda}')
else:
    print('  PyTorch CUDA: not available (CPU mode)')
" 2>/dev/null || echo "  PyTorch not available"
echo ""

# --- Validate inputs --------------------------------------------------------
if [[ ! -f "${ATAC_INPUT}" ]]; then
    echo "ERROR: SnapATAC2 AnnData not found: ${ATAC_INPUT}"
    echo "  Run script 02_snapatac2_processing.py first"
    exit 1
fi

if [[ ! -d "${RNA_DIR}" ]]; then
    echo "ERROR: CellRanger ARC directory not found: ${RNA_DIR}"
    echo "  Run script 01_cellranger_arc.sh first"
    exit 1
fi

# Create output directories
mkdir -p "${OUTPUT_DIR}/cistopic_model"
mkdir -p "${CHROMVAR_DIR}"
mkdir -p "${WORK_DIR}/logs"

# --- Phase 1: chromVAR TF motif analysis (Module 2c) -----------------------
echo ""
echo "============================================================"
echo "Phase 1: chromVAR TF Motif Analysis"
echo "  Started: $(date)"
echo "============================================================"

python "${SCRIPT_DIR}/03_chromvar_motifs.py" \
    --input "${ATAC_INPUT}" \
    --output-dir "${CHROMVAR_DIR}" \
    --genome hg38 \
    --n-background 50 \
    --n-top-heatmap 50

CHROMVAR_EXIT=$?
if [[ ${CHROMVAR_EXIT} -ne 0 ]]; then
    echo "WARNING: chromVAR analysis exited with code ${CHROMVAR_EXIT}"
    echo "  Continuing to SCENIC+ (chromVAR failure is non-fatal)"
fi

echo ""
echo "Phase 1 complete: $(date)"

# --- Phase 2: SCENIC+ GRN (Module 3) ---------------------------------------
echo ""
echo "============================================================"
echo "Phase 2: SCENIC+ Hepatocyte GRN"
echo "  Started: $(date)"
echo "============================================================"

python "${SCRIPT_DIR}/04_scenic_plus_grn.py" \
    --atac-input "${ATAC_INPUT}" \
    --rna-dir "${RNA_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-cells ${N_CELLS} \
    --n-topics "${N_TOPICS}" \
    --n-iter ${N_ITER} \
    --n-hvg ${N_HVG} \
    --genome hg38

SCENIC_EXIT=$?
if [[ ${SCENIC_EXIT} -ne 0 ]]; then
    echo "ERROR: SCENIC+ analysis exited with code ${SCENIC_EXIT}"
    exit ${SCENIC_EXIT}
fi

# --- Summary ----------------------------------------------------------------
echo ""
echo "============================================================"
echo "Pipeline COMPLETE"
echo "  Finished:   $(date)"
echo "  chromVAR:   ${CHROMVAR_DIR}/"
echo "  SCENIC+:    ${OUTPUT_DIR}/"
echo ""
echo "Output files:"
echo "  chromVAR:"
ls -lh "${CHROMVAR_DIR}/" 2>/dev/null || echo "    (no files)"
echo ""
echo "  SCENIC+:"
ls -lh "${OUTPUT_DIR}/"*.csv 2>/dev/null || echo "    (no files)"
echo "============================================================"
