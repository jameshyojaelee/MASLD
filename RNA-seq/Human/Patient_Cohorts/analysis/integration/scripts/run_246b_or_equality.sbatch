#!/bin/bash
#SBATCH --job-name=p25b_or_eq
#SBATCH --output=logs/p25b_or_eq_%j.log
#SBATCH --error=logs/p25b_or_eq_%j.log
#SBATCH --time=02:00:00
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -e
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Phase 2.5b: OR-equality bootstrap + look-elsewhere correction ==="
echo "Host: $(hostname)"
echo "Started: $(date)"

Rscript scripts/246b_or_equality_lookelsewhere.R

echo "Finished: $(date)"
