#!/bin/bash
#SBATCH --job-name=f3_substate_sc_donor
#SBATCH --output=logs/f3_substate_sc_donor_%j.log
#SBATCH --error=logs/f3_substate_sc_donor_%j.log
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

echo "=== F3 sub-state (sc) donor heterogeneity ==="
echo "Started: $(date)"

Rscript scripts/253_sc_donor_heterogeneity.R

echo "Finished: $(date)"
