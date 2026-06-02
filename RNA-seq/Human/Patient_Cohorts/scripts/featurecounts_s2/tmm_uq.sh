#!/bin/bash
#SBATCH --job-name=tmm_uq
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/logs/tmm_uq_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/logs/tmm_uq_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export R_PARALLEL_SEED=42

echo "=== TMM-vs-UQ normalization sensitivity (-s 2). Writes audit_sensitivity/tmm_vs_uq/dream_results_uq.csv ==="
echo "Started: $(date)  Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript RNA-seq/scripts/14_5_uq_normalization_sensitivity.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 14_5_uq_normalization_sensitivity.R"; exit 1; fi

echo "Finished: $(date)"
