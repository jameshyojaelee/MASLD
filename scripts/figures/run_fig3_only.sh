#!/bin/bash
#SBATCH --job-name=liver_fig3_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fig3_regen_%j.log

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs figures/main/fig2_genetics/panels

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Generating Fig 3: Causal Architecture ==="
Rscript scripts/figures/fig3_compact.R

echo ""
echo "=== Generating Supplementary Causal Extended ==="
Rscript scripts/figures/figS_causal_extended.R

echo ""
echo "Output:"
ls -la figures/main/fig2_genetics/genetics_compact.pdf 2>/dev/null || echo "fig3_compact.pdf not found"
ls -la figures/main/fig2_genetics/panels/*.pdf 2>/dev/null || echo "No individual panels"
ls -la figures/figS_causal_extended/*.pdf 2>/dev/null || echo "No supp panels"
