#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_obs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_obs_%j.err

# NOTE: no `set -u` — the rnaseq env's activate-binutils hook references unbound vars.
set -eo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$BASE"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript "$BASE/scripts/figures/fig2_eqtl_observability_boundary.R"
