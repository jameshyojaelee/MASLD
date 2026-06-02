#!/bin/bash
#SBATCH --job-name=sex_power
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_power_analysis/slurm_%j.log

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== 26b_sex_power_analysis.R ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"

micromamba run -n rnaseq Rscript \
  RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/26b_sex_power_analysis.R

echo "=== DONE: $(date) ==="
