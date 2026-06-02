#!/bin/bash
#SBATCH --job-name=composite_sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/57b_composite_%j.out
#SBATCH --error=logs/57b_composite_%j.err

set -e

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/GWAS/finemapping/src"

mkdir -p logs
mkdir -p "$PROJ/RNA-seq/results/audit_sensitivity/atac"

echo "[$(date)] Starting 57b_composite_score_sensitivity.R"
Rscript 57b_composite_score_sensitivity.R
echo "[$(date)] Done"
