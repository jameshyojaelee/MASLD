#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=gwas_atac_57_figures
#SBATCH --output=../logs/57_figures_%j.out
#SBATCH --error=../logs/57_figures_%j.err

source ~/.bashrc
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src"
cd "$SCRIPT_DIR"

echo "Running 57_gwas_atac_figures.R on $(hostname) at $(date)"
Rscript 57_gwas_atac_figures.R
echo "Done at $(date)"
