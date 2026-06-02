#!/bin/bash
#SBATCH --partition=bigmem --cpus-per-task=4 --mem=220G --time=12:00:00
#SBATCH --job-name=masld_figrefresh --output=RNA-seq/logs/masld_figrefresh_%j.out --error=RNA-seq/logs/masld_figrefresh_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== Script 209 (per-program COLOC + fig panels) ==="
Rscript RNA-seq/209_sex_subtype_figures.R 2>&1 | tail -40
echo "=== fig5 sex-dimorphism regen ==="
Rscript scripts/figures/fig5_sex_dimorphism.R 2>&1 | tail -30 || echo "fig5 failed (non-fatal)"
echo "=== fig3_progression regen ==="
Rscript scripts/figures/fig3_progression.R 2>&1 | tail -20 || echo "fig3 failed (non-fatal)"
echo "=== DONE ==="
