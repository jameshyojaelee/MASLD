#!/bin/bash
#SBATCH --job-name=fig_prognosis
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/fig_prognosis_%j.out
#SBATCH --error=scripts/figures/fig_prognosis_%j.err

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript scripts/figures/fig_prognosis_main.R

echo "=== DONE ==="
