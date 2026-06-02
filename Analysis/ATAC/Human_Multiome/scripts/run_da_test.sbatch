#!/bin/bash
#SBATCH --job-name=da_test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/da_test_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/da_test_%j.err

set -euo pipefail

ATAC_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome"

echo "DA Test - $(date)"
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

python "${ATAC_DIR}/scripts/07_run_da_test.py" \
    --input "${ATAC_DIR}/results/snapatac2/snapatac2_processed.h5ad" \
    --output-dir "${ATAC_DIR}/results/snapatac2"

echo "DA Test COMPLETE: $(date)"
