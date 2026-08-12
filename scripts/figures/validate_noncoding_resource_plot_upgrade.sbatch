#!/bin/bash
#SBATCH --job-name=ncplotval
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncplotval_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncplotval_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PYTHON="/gpfs/commons/home/jameslee/micromamba/envs/alphagenome/bin/python"
SOURCE="${PROJECT_ROOT}/figures/candidates/noncoding-resource-2026-08-12/plot-upgrade-v3"
REVIEW="${PROJECT_ROOT}/figures/candidates/noncoding-resource-2026-08-12/plot-upgrade-v3-review"

test -d "${SOURCE}"
test ! -e "${REVIEW}"
"${PYTHON}" scripts/figures/validate_noncoding_resource_plot_upgrade.py \
  --source "${SOURCE}" \
  --review "${REVIEW}"
