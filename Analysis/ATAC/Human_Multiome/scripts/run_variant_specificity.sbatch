#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_Multiome/logs/variant_specificity_%j.out
#SBATCH --error=Analysis/ATAC/Human_Multiome/logs/variant_specificity_%j.err
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
eval "$(micromamba shell hook -s bash)"; micromamba activate snapatac2
python Analysis/ATAC/Human_Multiome/scripts/58_variant_celltype_specificity.py
echo "=== Done $(date) ==="
