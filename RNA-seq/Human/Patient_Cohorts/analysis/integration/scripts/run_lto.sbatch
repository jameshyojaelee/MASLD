#!/bin/bash
#SBATCH --job-name=LTO
#SBATCH --partition=cpu
#SBATCH --array=1-10%5
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=logs/LTO_%A_%a.out
#SBATCH --error=logs/LTO_%A_%a.err

# Pillar A — Leave-TWO-cohorts-out (10 pairs from 5 mega cohorts).
# Cohort-level stress test at ~60% retention. Complements LOCO at 80%.

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export LTO_PAIR_IDX=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] LTO pair ${LTO_PAIR_IDX} on $(hostname)"
Rscript dream_lto_iter.R
echo "[$(date)] Done LTO pair ${LTO_PAIR_IDX}"
