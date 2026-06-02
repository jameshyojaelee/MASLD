#!/bin/bash
#SBATCH --job-name=14p4_arr
#SBATCH --array=1-150
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/14p4_arr_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/14p4_arr_%A_%a.err

# 14.4_sex_subsampling_iter array
# Map array task id 1-150 to (FRAC, ITER):
#   FRAC = c(0.50, 0.60, 0.70, 0.80, 0.90)[((id-1) %/% 30) + 1]
#   ITER = ((id-1) %% 30) + 1
#
# Per-task env vars: SCRIPT controls which iter script runs (14.4 vs 14.4c).
# Submit two arrays (one per SCRIPT) using --export.

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

ID="${SLURM_ARRAY_TASK_ID}"
FRACS=(0.50 0.60 0.70 0.80 0.90)
FRAC_IDX=$(( (ID - 1) / 30 ))
ITER=$(( (ID - 1) % 30 + 1 ))
export SUBSAMPLE_FRAC=${FRACS[$FRAC_IDX]}
export SUBSAMPLE_ITER=$ITER

SCRIPT_NAME=${SCRIPT:-14.4_sex_subsampling_iter.R}

echo "Started: $(date)  TASK=${ID}  FRAC=${SUBSAMPLE_FRAC}  ITER=${SUBSAMPLE_ITER}  SCRIPT=${SCRIPT_NAME}  CPUs=${SLURM_CPUS_PER_TASK}"
Rscript "${SCRIPT_NAME}"
echo "Finished: $(date)"
