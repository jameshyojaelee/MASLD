#!/bin/bash
#SBATCH --job-name=sc_cohort_loco_cv
#SBATCH --output=logs/sc_cohort_loco_cv_%j.log
#SBATCH --error=logs/sc_cohort_loco_cv_%j.log
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -e
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== sc cohort LOCO-CV bimodality (253c) ==="
echo "Started: $(date)"

Rscript scripts/253c_sc_cohort_loco_cv.R

echo "Finished: $(date)"
