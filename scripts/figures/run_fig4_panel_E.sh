#!/bin/bash
#SBATCH --job-name=fig4_E_drug
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig4_E_drug_%j.out
#SBATCH --error=logs/fig4_E_drug_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/main/fig4_validation/panels
Rscript scripts/figures/fig4_panel_E_drug_regulon_disruption.R
