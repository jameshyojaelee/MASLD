#!/bin/bash
#SBATCH --job-name=liver_fig1_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fig1_regen_%j.log

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs figures

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Generate panel versions of embedded PDFs (standalone originals unchanged)
python scripts/figures/fig_sunburst.py
python scripts/figures/fig_convergence_wheel.py --no-title
python scripts/figures/fig_bulkrna_matrix.py --panel

Rscript scripts/figures/fig1_compact.R
ls -la figures/fig1/fig1_compact.pdf
