#!/bin/bash
#SBATCH --job-name=figS1_volcano
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/figS1_volcano_%j.out
#SBATCH --error=scripts/figures/logs/figS1_volcano_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "[$(date)] Starting figS_deg_validation.R..."
Rscript scripts/figures/figS_deg_validation.R
echo "[$(date)] Done."

echo ""
echo "Output files:"
ls -lh figures/supplementary/figS02_progression/figS_deg_validation.pdf 2>/dev/null
echo ""
echo "Panel files:"
ls -lh figures/supplementary/figS02_progression/ 2>/dev/null
