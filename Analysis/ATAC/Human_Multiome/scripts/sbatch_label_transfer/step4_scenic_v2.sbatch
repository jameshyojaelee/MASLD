#!/bin/bash
#SBATCH --job-name=scenic_v2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step4_scenic_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step4_scenic_v2_%j.err

set -euo pipefail
export ATAC_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome"
export SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/scripts"
export PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export RNA_ATLAS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
export DONOR_PAIRING="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/GSE244832/metadata/donor_pairing.csv"
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate atac_env
CONDA_PREFIX_PATH="$(micromamba info | grep 'env location' | awk '{print $NF}')"
if [[ -d "${CONDA_PREFIX_PATH}/lib" ]]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX_PATH}/lib:${LD_LIBRARY_PATH:-}"
fi

cd "$ATAC_DIR"
RELABELED="results/label_transfer/snapatac2_relabeled.h5ad"
if [[ ! -f "${RELABELED}" ]]; then
    echo "ERROR: relabeled h5ad not found. Step 2 may have failed."
    exit 1
fi

python "$SCRIPT_DIR/04_scenic_plus_grn.py" \
    --atac-input "${RELABELED}" \
    --rna-dir cellranger_arc \
    --output-dir results/scenic_plus_v2 \
    --n-cells 15000 \
    --n-topics "10,20,30,40" \
    --n-iter 500 \
    --n-hvg 3000 \
    --genome hg38

echo "Step 4 COMPLETE: SCENIC+ re-run done"
