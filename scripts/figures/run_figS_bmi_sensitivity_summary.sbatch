#!/bin/bash
#SBATCH --job-name=figS_bmi_sens
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS_bmi_sensitivity_summary_%j.log

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript scripts/figures/figS_bmi_sensitivity_summary.R
