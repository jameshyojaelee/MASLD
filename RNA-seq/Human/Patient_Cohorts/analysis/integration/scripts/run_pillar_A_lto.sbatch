#!/bin/bash
#SBATCH --job-name=pillarA_lto
#SBATCH --partition=cpu
#SBATCH --array=1-10%5
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarA_lto_%A_%a.out
#SBATCH --error=logs/pillarA_lto_%A_%a.err

# Pillar A — Leave-TWO-cohorts-out (10 pairs from 5 mega cohorts).
# Cohort-level stress test at ~60 % retention. Complements LOCO at 80 %.

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export LTO_PAIR_IDX=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] Starting LTO pair ${LTO_PAIR_IDX} on $(hostname)"
Rscript dream_lto_iter.R
echo "[$(date)] Done LTO pair ${LTO_PAIR_IDX}"
