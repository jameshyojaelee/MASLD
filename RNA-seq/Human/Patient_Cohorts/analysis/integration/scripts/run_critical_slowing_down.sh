#!/bin/bash
#SBATCH --job-name=critical_slowing_down
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/critical_slowing_down_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/critical_slowing_down_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs

echo "=== Critical Slowing Down Analysis ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript 195_critical_slowing_down.R

echo "End: $(date)"
