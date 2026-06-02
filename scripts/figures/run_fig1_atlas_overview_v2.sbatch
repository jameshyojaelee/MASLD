#!/bin/bash
#SBATCH --job-name=fig1_atlas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig1_atlas_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig1_atlas_v2_%j.err

set -eo pipefail
source "$HOME/.bashrc"
micromamba activate rnaseq

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$BASE/scripts/figures/logs"

cd "$BASE"
Rscript "$BASE/scripts/figures/fig1_atlas_overview_v2.R"
