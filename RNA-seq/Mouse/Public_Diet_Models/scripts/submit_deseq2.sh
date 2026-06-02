#!/bin/bash
#SBATCH --job-name=deseq2_varied
#SBATCH --output=logs/deseq2_%j.out
#SBATCH --error=logs/deseq2_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

# Run DESeq2 analysis for all varied diet mouse datasets
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq

# Use unified_env which has all DESeq2 dependencies
MAMBA_ENV="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"

echo "Starting DESeq2 analysis at $(date)"
micromamba run -p $MAMBA_ENV Rscript scripts/run_deseq2_varied_datasets.R
echo "Completed DESeq2 analysis at $(date)"
