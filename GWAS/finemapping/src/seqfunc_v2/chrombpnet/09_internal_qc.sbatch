#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-4%5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/internal_qc_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/internal_qc_%A_%a.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python3 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet/09_internal_qc.py" \
  --fold "${SLURM_ARRAY_TASK_ID}" \
  --require-pass
