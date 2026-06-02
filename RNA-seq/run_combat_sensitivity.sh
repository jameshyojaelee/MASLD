#!/bin/bash
#SBATCH --job-name=liver_combat_sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=logs/combat_sensitivity_%j.log

set -euo pipefail

echo "=== ComBat-seq Sensitivity Analysis ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 34_combat_seq_sensitivity.R

echo "End: $(date)"
