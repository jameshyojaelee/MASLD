#!/bin/bash
#SBATCH --job-name=mpra
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/56j_zhu_%j.out
#SBATCH --error=logs/56j_zhu_%j.err

set -eo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
mkdir -p GWAS/finemapping/src/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[$(date)] Starting Zhu/Hu 2026 MPRA/CRISPRi concordance (job $SLURM_JOB_ID)"
echo "Host: $(hostname)"

Rscript GWAS/finemapping/src/56j_zhu_mpra_concordance.R

echo "[$(date)] Done."
