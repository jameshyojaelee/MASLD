#!/bin/bash -l
#SBATCH --job-name=topld_test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-22
#SBATCH --output=GWAS/finemapping/logs/topld_test_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/topld_test_%A_%a.err

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

CHR=${SLURM_ARRAY_TASK_ID}
export LD_PANEL=topld
export COLOC_OUT_SUFFIX="_topld"

echo "[topld_test] UKBB_GGT chr${CHR} LD_PANEL=${LD_PANEL} start $(date)"
Rscript GWAS/finemapping/src/06_susie_coloc.R UKBB_GGT ${CHR}
echo "[topld_test] UKBB_GGT chr${CHR} done $(date)"
