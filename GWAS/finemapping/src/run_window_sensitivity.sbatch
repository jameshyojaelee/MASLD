#!/bin/bash
#SBATCH --job-name=coloc_window_sens
#SBATCH --output=logs/coloc_window_sensitivity_%j.out
#SBATCH --error=logs/coloc_window_sensitivity_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=48:00:00

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

echo "Job: $SLURM_JOB_ID | Node: $SLURM_NODELIST | Start: $(date)"

micromamba run -n rnaseq Rscript src/06b_coloc_window_sensitivity.R

echo "Done: $(date)"
