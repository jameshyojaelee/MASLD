#!/bin/bash
#SBATCH --job-name=figS_sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=scripts/figures/logs/figS_sensitivity_%j.out
#SBATCH --error=scripts/figures/logs/figS_sensitivity_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript scripts/figures/figS_sensitivity.R

echo "=== Done: $(date) ==="
