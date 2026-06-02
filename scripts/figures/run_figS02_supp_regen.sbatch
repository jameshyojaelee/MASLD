#!/bin/bash
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --job-name=figS02_supp_regen
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figS02_supp_regen_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figS02_supp_regen_%j.err
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

set -e
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

Rscript scripts/figures/figS_disease_progression_supp.R || echo "figS_disease_progression_supp.R failed"
Rscript scripts/figures/figS_progression_landscape.R || echo "figS_progression_landscape.R failed"
