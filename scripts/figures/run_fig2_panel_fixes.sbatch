#!/bin/bash
#SBATCH --job-name=fig2_fix
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_panel_fixes_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_panel_fixes_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

module purge 2>/dev/null || true
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Regenerating Fig 2 (12 panels with fixes) ==="
echo "Start: $(date)"

Rscript scripts/figures/fig2_compact.R

echo "=== Done: $(date) ==="
echo "--- Checking outputs ---"
ls -la figures/fig2/fig2_compact.pdf
ls -la figures/fig2/panels/fig2*_*.pdf 2>/dev/null || echo "No panel PDFs found"
