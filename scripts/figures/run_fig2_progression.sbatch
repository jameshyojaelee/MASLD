#!/bin/bash
#SBATCH --job-name=fig2_progression
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/fig2_progression_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/fig2_progression_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting Fig 2 panel renders + compositor on $(hostname)"

# Render each new panel
for s in fig2_panel_cascade_degs fig2_panel_inflection_composite \
         fig2_panel_celltype_cascade fig2_panel_network_communities \
         fig2_panel_nmf_programs; do
  echo "[$(date)] panel script: $s"
  micromamba run -n rnaseq Rscript scripts/figures/$s.R
done

# Compose
echo "[$(date)] compositor"
micromamba run -n rnaseq Rscript scripts/figures/fig2_progression.R

echo "[$(date)] DONE"
