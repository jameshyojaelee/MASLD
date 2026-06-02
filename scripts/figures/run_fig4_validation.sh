#!/bin/bash
#SBATCH --job-name=fig4_val
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/fig4_val_%j.out
#SBATCH --error=logs/fig4_val_%j.err

set -eo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$BASE/scripts/figures"
mkdir -p logs

# shellcheck disable=SC1091
source "$HOME/.bashrc"
micromamba activate rnaseq
set -u

Rscript fig4_validation.R
