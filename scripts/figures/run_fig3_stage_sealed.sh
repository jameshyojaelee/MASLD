#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3_stage_sealed_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3_stage_sealed_%j.log

# Renders the cross-sectional adjacent-stage Fig 3 panels from the sealed
# fibrosis-adjacent-true-kleiner-lvqw-v1 bundle. The R script verifies the
# bundle SHA256 and asserts every count before drawing anything.

# NB: no `set -u` — it breaks the rnaseq env activation hook
# (activate-binutils_linux-64.sh references unbound ADDR2LINE).
set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "${BASE}/scripts/figures/logs"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript "${BASE}/scripts/figures/fig3_stage_degs_sealed_crosssectional.R"
