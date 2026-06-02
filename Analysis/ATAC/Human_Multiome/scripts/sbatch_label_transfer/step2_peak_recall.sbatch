#!/bin/bash
#SBATCH --job-name=peak_recall
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step2_peak_recall_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step2_peak_recall_%j.err

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
python "$SCRIPT_DIR/02c_rerun_peaks_with_transferred_labels.py" \
    --label_transferred results/label_transfer/snapatac2_label_transferred.h5ad \
    --per_donor_dir results/snapatac2/per_donor \
    --old_peak_dir results/snapatac2/cell_type_peak_sets \
    --output_dir results/label_transfer

echo "Step 2 COMPLETE: Peak re-calling done"
