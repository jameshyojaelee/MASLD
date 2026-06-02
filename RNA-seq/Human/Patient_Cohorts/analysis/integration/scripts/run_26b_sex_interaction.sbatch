#!/bin/bash
#SBATCH --job-name=m2_26b_sex_int
#SBATCH --output=logs/26b_sex_int_%j.log
#SBATCH --error=logs/26b_sex_int_%j.log
#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

set -e
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Team M2: 26b sex x stage x transition interaction ==="
echo "Started: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo

Rscript scripts/26b_sex_interaction_per_transition.R

echo
echo "=== 26b complete: $(date) ==="
