#!/bin/bash
#SBATCH --job-name=robust_agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --output=logs/robust_agg_%j.out
#SBATCH --error=logs/robust_agg_%j.err

# Aggregate Pillar A / C / D and assemble the master robustness atlas.
# Submit with afterok dependencies on all the array jobs:
#   sbatch --dependency=afterok:${CPSS_JOB}:${BOOT_JOB}:${KFOLD_JOB}:${PERM_JOB}:${PVCA_JOB}:${SVA_JOB}:${PILLAR_B_JOB} run_aggregate_robustness.sbatch

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Starting aggregation on $(hostname)"

Rscript aggregate_pillar_A.R
Rscript aggregate_pillar_C.R
Rscript aggregate_robustness_atlas.R

echo "[$(date)] Aggregation complete"
