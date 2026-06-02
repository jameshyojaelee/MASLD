#!/bin/bash
#SBATCH --job-name=format_bbj
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/BBJ/logs/format_bbj_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/BBJ/logs/format_bbj_%j.err

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "Formatting BBJ GWAS for COLOC (hg19 → hg38 liftover)"
echo "Started: $(date)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/BBJ
Rscript format_bbj_for_coloc.R 2>&1

echo "Completed: $(date)"
