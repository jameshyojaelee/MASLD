#!/bin/bash
#SBATCH --job-name=127_spatial_val
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/127_spatial_validation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/127_spatial_validation_%j.err

set -euo pipefail

echo "=== Script 127: Spatial Validation of Progression Signatures ==="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

SCRIPT_DIR=${MASLD_PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

python3 ${SCRIPT_DIR}/127_spatial_validation.py

echo "=== Completed: $(date) ==="
