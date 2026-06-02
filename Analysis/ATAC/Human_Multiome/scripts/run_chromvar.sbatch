#!/bin/bash
#SBATCH --job-name=chromVAR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/chromvar_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/chromvar_%j.err

# =============================================================================
# Module 2, Step 3: chromVAR TF motif enrichment
# Depends on SnapATAC2 processing (02_snapatac2_processing.py)
# =============================================================================

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATAC_DIR="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome"

echo "============================================================"
echo "chromVAR TF Motif Analysis"
echo "  Job ID:     ${SLURM_JOB_ID}"
echo "  Node:       $(hostname)"
echo "  Started:    $(date)"
echo "============================================================"

module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# Use fixed h5ad (unique obs_names) if available, else original
INPUT="${ATAC_DIR}/results/snapatac2/snapatac2_processed_fixed.h5ad"
if [[ ! -f "$INPUT" ]]; then
    INPUT="${ATAC_DIR}/results/snapatac2/snapatac2_processed.h5ad"
fi
echo "  Input: $INPUT"

python "${ATAC_DIR}/scripts/03_chromvar_motifs.py" \
    --input "$INPUT" \
    --output-dir "${ATAC_DIR}/results/chromvar" \
    --genome hg38

echo ""
echo "chromVAR COMPLETE: $(date)"
