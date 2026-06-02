#!/bin/bash
#SBATCH --job-name=fig1g_lfc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig1g_lfc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig1g_lfc_%j.err

set -eo pipefail
source "$HOME/.bashrc"
micromamba activate rnaseq

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$BASE/scripts/figures/logs"

cd "$BASE"
Rscript scripts/figures/fig1g_patient_lfc_combined.R
