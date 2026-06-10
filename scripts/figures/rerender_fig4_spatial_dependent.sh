#!/usr/bin/env bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/rerender_fig4_spatial_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/rerender_fig4_spatial_%j.out

# Re-render the figures that depend on the spatial-integration atlas (rebuilt by 06)
# + fig4_validation (Panel 4c SVG-schema fix). NOTE: do NOT use `set -u` — micromamba's
# binutils activation hook references an unbound ADDR2LINE and trips it (matches run_pub_figures.sh).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures
mkdir -p logs
echo "host=$(hostname) start=$(date)"
for s in fig4a_proteomics_overview.R fig4f_fads2_zonation.R fig4f_cyp3a4_zonation.R fig4_validation.R; do
  echo "=== render $s ==="
  Rscript "$s" 2>&1 || echo "WARN: $s rc=$?"
done
echo "=== DONE $(date) ==="
