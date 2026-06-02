#!/bin/bash
#SBATCH --job-name=figS13_ccc
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figS13_ccc_dynamics_%j.log

set -euo pipefail

source ~/.bashrc
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

Rscript "${MASLD_PROJECT_ROOT}/scripts/figures/figS_ccc_dynamics.R"
