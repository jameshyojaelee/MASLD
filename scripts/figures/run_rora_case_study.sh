#!/bin/bash
#SBATCH --job-name=rora_case_study
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/rora_case_study_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/rora_case_study_%j.err

set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

mkdir -p "${BASE}/scripts/figures/logs"

echo "[rora_case_study] start $(date)"
cd "${BASE}"
Rscript scripts/figures/fig_rora_case_study.R
echo "[rora_case_study] done  $(date)"
