#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_resolution_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_resolution_%j.err

# NB: no `set -u` -- the rnaseq activate.d hook references unbound vars
# (ADDR2LINE) and aborts under nounset.
set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript "${BASE}/scripts/figures/fig4_convergence_resolution_limits.R"
