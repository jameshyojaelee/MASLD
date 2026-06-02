#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=1-22
#SBATCH --job-name=liver-coloc
#SBATCH --output=logs/coloc_cw2_%A_%a.out
#SBATCH --error=logs/coloc_cw2_%A_%a.err

# SuSiE-COLOC for coworker2 — uses full env path (no local conda needed)
# Usage: GWAS_NAME=Ghouse_HCC sbatch src/06_susie_coloc_coworker2.sh

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/finemapping/bin/Rscript"

export CHR_FILTER="${SLURM_ARRAY_TASK_ID}"

echo "=== SuSiE-COLOC (coworker2): ${GWAS_NAME} chr${CHR_FILTER} ==="
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "User: $(whoami)"
echo "Start: $(date)"

${RSCRIPT} src/06_susie_coloc.R "${GWAS_NAME}" "${CHR_FILTER}"

# Cleanup shared checkpoint
rm -f "results/susie_coloc/${GWAS_NAME}/checkpoint_chr${CHR_FILTER}.csv" 2>/dev/null || true

echo "End: $(date)"
