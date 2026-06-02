#!/bin/bash
#SBATCH --job-name=fig2_rev
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/fig2_rev_%j.log

set -eo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript scripts/figures/fig2_progression_sex.R

echo "=== Outputs ==="
ls -la figures/main/fig2_progression_sex/ figures/main/fig2_progression_sex/panels/
