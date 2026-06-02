#!/bin/bash
#SBATCH --job-name=currin_caqtl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/56i_currin_%j.out
#SBATCH --error=logs/56i_currin_%j.err

set -eo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
mkdir -p GWAS/finemapping/src/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[$(date)] Starting Currin caQTL B4 concordance (job $SLURM_JOB_ID)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

Rscript GWAS/finemapping/src/56i_currin_caqtl_concordance.R

echo "[$(date)] Done."
