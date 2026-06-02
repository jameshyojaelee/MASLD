#!/bin/bash
#SBATCH --job-name=label_transfer
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step1_label_transfer_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step1_label_transfer_%j.err

set -euo pipefail
export ATAC_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome"
export SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/scripts"
export PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export RNA_ATLAS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
export DONOR_PAIRING="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/GSE244832/metadata/donor_pairing.csv"
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd "$ATAC_DIR"
python "$SCRIPT_DIR/02b_label_transfer_from_rna.py" \
    --atac_dir results/snapatac2 \
    --rna_atlas "$RNA_ATLAS" \
    --donor_pairing "$DONOR_PAIRING" \
    --output_dir results/label_transfer \
    --k_neighbors 10 \
    --n_hvgs 3000 \
    --n_pcs 30

echo "Step 1 COMPLETE: Label transfer done"
