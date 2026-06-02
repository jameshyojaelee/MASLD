#!/bin/bash
#SBATCH --job-name=A7_T3_perm
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=90:00:00
#SBATCH --array=1-10
#SBATCH --output=logs/A7_T3_perm_%A_%a.out
#SBATCH --error=logs/A7_T3_perm_%A_%a.err

set -eo pipefail
mkdir -p logs

echo "=== A7 T3: Cohort-stratified sex permutation null ==="
echo "Job ID: ${SLURM_ARRAY_JOB_ID:-NA}  Task: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Node: $(hostname)"
echo "Started: $(date)"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export R_PARALLEL_SEED=42
export N_PERMS_PER_CHUNK=50
export B_TOTAL=500

Rscript 14.5a_T3_cohort_perm.R

echo "Task ${SLURM_ARRAY_TASK_ID:-NA} done: $(date)"
