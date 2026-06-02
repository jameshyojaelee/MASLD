#!/bin/bash
#SBATCH --job-name=CPSS
#SBATCH --partition=cpu
#SBATCH --array=1-100%25
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=logs/CPSS_%A_%a.out
#SBATCH --error=logs/CPSS_%A_%a.err

# Pillar A — Complementary-Pairs Stability Selection (B=100 pairs => 200 fits).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export CPSS_PAIR=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] CPSS pair ${CPSS_PAIR} on $(hostname)"
Rscript dream_cpss_iter.R
echo "[$(date)] Done CPSS pair ${CPSS_PAIR}"
