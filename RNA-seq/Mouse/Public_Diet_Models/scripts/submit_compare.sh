#!/bin/bash
#SBATCH --job-name=compare_datasets
#SBATCH --output=logs/compare_datasets_%j.out
#SBATCH --error=logs/compare_datasets_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq

MAMBA_ENV="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"

echo "Starting comparative analysis at $(date)"
micromamba run -p $MAMBA_ENV Rscript scripts/compare_diet_vs_mcd.R
echo "Completed at $(date)"
