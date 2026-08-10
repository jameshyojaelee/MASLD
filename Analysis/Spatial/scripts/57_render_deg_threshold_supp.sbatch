#!/usr/bin/env bash
#SBATCH --job-name=figures
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/deg_threshold_supp_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/deg_threshold_supp_%j.out
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n spatial python scripts/figures/figS_deg_threshold_sensitivity_data.py
micromamba run -n rnaseq  Rscript scripts/figures/figS_deg_threshold_sensitivity.R
echo "EXIT=$?"
