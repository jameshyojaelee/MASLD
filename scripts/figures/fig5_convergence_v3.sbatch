#!/bin/bash
#SBATCH --job-name=fig5_v3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig5_v3_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig5_v3_%j.err

set -eo pipefail

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs

# shellcheck disable=SC1091
source "$HOME/.bashrc"
set +u
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript scripts/figures/fig5_convergence_v3.R
