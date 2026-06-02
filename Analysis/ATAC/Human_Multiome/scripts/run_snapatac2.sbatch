#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=36:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/snapatac2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/snapatac2_%j.err

# =============================================================================
# SnapATAC2 scATAC-seq Processing — GSE244832 Human Liver
# Module 2, Step 2: Fragment QC → Tile → LSI → Harmony → Clustering → DA
#
# Processes fragment files (18 donors) from custom alignment pipeline.
# GSE244832 uses combinatorial indexing (NOT 10x), so fragments are generated
# by 01_align_scatac.sh (bowtie2 + barcode extraction).
#
# Prerequisites:
#   1. Alignment completed (01_align_scatac.sh)
#   2. snapatac2 conda environment available
#
# Usage:
#   sbatch scripts/run_snapatac2.sbatch
# =============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATAC_DIR="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome"
SCRIPT="${ATAC_DIR}/scripts/02_snapatac2_processing.py"
INPUT_DIR="${ATAC_DIR}/results/fragments"
OUTPUT_DIR="${ATAC_DIR}/results/snapatac2"
SAMPLE_CSV="${PROJECT_ROOT}/data/GSE244832/metadata/donor_pairing.csv"

echo "============================================================"
echo "SnapATAC2 scATAC-seq Processing"
echo "  Job ID:     ${SLURM_JOB_ID}"
echo "  Node:       $(hostname)"
echo "  CPUs:       ${SLURM_CPUS_PER_TASK}"
echo "  Memory:     ${SLURM_MEM_PER_NODE}MB"
echo "  Started:    $(date)"
echo "============================================================"

# Activate environment with isolation fixes
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# Verify key packages
python -c "import snapatac2 as snap; print(f'SnapATAC2 version: {snap.__version__}')"
python -c "import anndata; print(f'AnnData version: {anndata.__version__}')"
macs3 --version

# --- Create output directories -----------------------------------------------
mkdir -p "${OUTPUT_DIR}"

# --- Run pipeline ------------------------------------------------------------
echo ""
echo "Running SnapATAC2 pipeline..."
echo "  Input:      ${INPUT_DIR}"
echo "  Output:     ${OUTPUT_DIR}"
echo "  Sample CSV: ${SAMPLE_CSV}"
echo ""

python "${SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --sample_csv "${SAMPLE_CSV}" \
    --n_workers "${SLURM_CPUS_PER_TASK}" \
    --min_tsse 4.0 \
    --min_frags 1000 \
    --max_frags 100000

echo ""
echo "============================================================"
echo "SnapATAC2 processing COMPLETE"
echo "  Finished:   $(date)"
echo "  Output:     ${OUTPUT_DIR}"
echo "============================================================"
