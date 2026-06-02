#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=gwas_atac_04_aggregate
#SBATCH --output=../logs/04_aggregate_%j.out
#SBATCH --error=../logs/04_aggregate_%j.err

source ~/.bashrc
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src"
cd "$SCRIPT_DIR"

echo "Running 04_aggregate_results.R on $(hostname) at $(date)"
Rscript 04_aggregate_results.R
echo "Done at $(date)"
