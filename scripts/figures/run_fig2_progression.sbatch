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

echo "[$(date)] starting Fig 2 panel renders on $(hostname)"

# Render each new panel
for s in cascade_degs inflection_composite \
         celltype_cascade network_communities \
         nmf_programs; do
  echo "[$(date)] panel script: $s"
  micromamba run -n rnaseq Rscript scripts/figures/$s.R
done

echo "[$(date)] DONE"
