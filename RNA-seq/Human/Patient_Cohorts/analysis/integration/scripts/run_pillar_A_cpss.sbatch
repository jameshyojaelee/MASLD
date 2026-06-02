#!/bin/bash
#SBATCH --job-name=pillarA_cpss
#SBATCH --partition=cpu
#SBATCH --array=1-100%30
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarA_cpss_%A_%a.out
#SBATCH --error=logs/pillarA_cpss_%A_%a.err

# Pillar A — Complementary-Pairs Stability Selection (B=100 pairs => 200 fits).
# Each task fits dream on a 50% half-sample I AND its complement I^c, with
# cohort-stratified subsampling. Throttled to 30 concurrent jobs to stay under
# the nslab 7TB QOS cap (30 * ~100GB = 3TB peak).

set -eo pipefail
mkdir -p logs

source ~/.bashrc 2>/dev/null || true
micromamba activate rnaseq

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
export CPSS_PAIR=${SLURM_ARRAY_TASK_ID}
echo "[$(date)] Starting CPSS pair ${CPSS_PAIR} on $(hostname)"
Rscript dream_cpss_iter.R
echo "[$(date)] Done CPSS pair ${CPSS_PAIR}"
