#!/bin/bash
#SBATCH --job-name=steatosis_bmi_sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/steatosis_bmi_sensitivity_%j.log
#SBATCH --error=logs/steatosis_bmi_sensitivity_%j.log

# NOTE: no `-u` (nounset) — micromamba/conda activate.d scripts (binutils ADDR2LINE)
# reference unset vars and abort the job under set -u before Rscript runs.
set -eo pipefail

echo "=== Steatosis/BMI Sensitivity Analysis ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript scripts/figures/figS_steatosis_bmi_sensitivity.R

echo "=== Completed: $(date) ==="
