#!/bin/bash
#SBATCH --job-name=decoupleR
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/51b_decoupler_sensitivity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/51b_decoupler_sensitivity_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# 51b_decoupler_collectri_sensitivity: DoRothEA vs CollecTRI TF activity comparison
# Addresses reviewer R2 section 6.C
# Submit: sbatch run_51b_decoupler_sensitivity.sh

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

echo "=== 51b: decoupleR DoRothEA vs CollecTRI Sensitivity ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo ""

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 51b_decoupler_collectri_sensitivity.R 2>&1

echo ""
echo "=== Completed: $(date) ==="
