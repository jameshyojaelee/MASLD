#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=gwas_atac_55_overlap
#SBATCH --output=../logs/55_overlap_%j.out
#SBATCH --error=../logs/55_overlap_%j.err

source ~/.bashrc
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src"
cd "$SCRIPT_DIR"

echo "Running 55_gwas_atac_variant_overlap.R on $(hostname) at $(date)"
Rscript 55_gwas_atac_variant_overlap.R
echo "Done at $(date)"
