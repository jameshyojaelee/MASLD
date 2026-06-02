#!/usr/bin/env bash
#SBATCH --job-name=fig4_spatial
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_spatial_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_spatial_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}"

# Main Figure 4 spatial panels (f)-(h)
Rscript scripts/figures/fig4_spatial_panels.R

echo "Done: fig4_spatial_panels.pdf + individual panels in figures/fig4_spatial_panels/"
