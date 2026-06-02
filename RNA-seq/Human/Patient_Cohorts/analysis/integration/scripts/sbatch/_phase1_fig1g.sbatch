#!/bin/bash
#SBATCH --job-name=p1_fig1g
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p1_fig1g_%j.out
#SBATCH --error=logs/p1_fig1g_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Started: $(date)"
Rscript scripts/figures/fig1g_loo_lfc_stability.R
echo "Finished: $(date)"
