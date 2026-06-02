#!/bin/bash
#SBATCH --job-name=p1_6_power_analysis
#SBATCH --output=logs/p1_6_power_analysis_%j.log
#SBATCH --error=logs/p1_6_power_analysis_%j.log
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -e
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Phase 1.6 Power Analysis ==="
echo "Host: $(hostname)"
echo "Started: $(date)"

Rscript scripts/240_power_analysis.R

echo "Finished: $(date)"
