#!/bin/bash
#SBATCH --job-name=fig2_figs_only
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2_figs_only_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2_figs_only_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"
cd "${BASE}"

echo "=== Figure 2 + FigS8 Generation: $(date) ==="

echo "=== Figure 2: Disease Progression ==="
Rscript scripts/figures/fig2_disease_progression.R 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: fig2_disease_progression.R"
  exit 1
fi

echo ""
echo "=== Supplementary Figure 8 ==="
Rscript scripts/figures/figS8_disease_progression_supp.R 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: figS8_disease_progression_supp.R"
  exit 1
fi

echo ""
echo "=== All figures completed: $(date) ==="
ls -lh "${BASE}/figures/fig2_disease_progression.pdf" \
       "${BASE}/figures/figS8_disease_progression_supp.pdf" \
       2>/dev/null || true
