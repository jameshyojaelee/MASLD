#!/bin/bash
#SBATCH --job-name=fig1d_ideas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=figures/logs/fig1d_ideas_%j.out
#SBATCH --error=figures/logs/fig1d_ideas_%j.err

set -euo pipefail
BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
export MASLD_PROJECT_ROOT="${BASE}"
cd "${BASE}"
mkdir -p figures/logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript "${BASE}/scripts/figures/fig1d_candidates.R"
