#!/bin/bash
#SBATCH --job-name=sex_int_sub
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --array=1-30
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling/logs/sub_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling/logs/sub_%A_%a.err

# ---------------------------------------------------------------------------
# SLURM array: 30 iterations of interaction-based sex subsampling
# Each array task runs ONE iteration at 70% subsampling fraction
#
# Usage:
#   sbatch run_sex_interaction_subsampling.sh
#   # After all array tasks complete, run aggregation:
#   sbatch --dependency=afterok:$SLURM_ARRAY_JOB_ID run_sex_interaction_aggregate.sh
# ---------------------------------------------------------------------------

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.4c_sex_interaction_subsampling_iter.R"
LOG_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling/logs"

mkdir -p "$LOG_DIR"

export SUBSAMPLE_FRAC=0.70
export SUBSAMPLE_ITER=${SLURM_ARRAY_TASK_ID}

echo "============================================================"
echo "  Sex Interaction Subsampling: FRAC=${SUBSAMPLE_FRAC} ITER=${SUBSAMPLE_ITER}"
echo "  Date: $(date)"
echo "  Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Node: $(hostname)"
echo "  CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"

Rscript "$SCRIPT"

echo "============================================================"
echo "  DONE: $(date)"
echo "============================================================"
