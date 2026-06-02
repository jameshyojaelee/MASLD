#!/bin/bash
#SBATCH --job-name=sex_consolidation_regen
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_consolidation_regen_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_consolidation_regen_%j.err

set -eo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p "${BASE}/figures/logs"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "============================================================"
echo "  Sex stratified-causal consolidation regen"
echo "  Date: $(date)"
echo "  Job ID: ${SLURM_JOB_ID:-local}"
echo "============================================================"

echo ""
echo "--- 1. figS_coloc_sex_direct.R (csd_a/b/c/e/f + composite) ---"
Rscript "${BASE}/scripts/figures/figS_coloc_sex_direct.R" || \
  echo "WARN: figS_coloc_sex_direct.R failed (rc=$?)"

echo ""
echo "--- 2. RNA-seq/209_sex_subtype_figures.R ---"
Rscript "${BASE}/RNA-seq/209_sex_subtype_figures.R" || \
  echo "WARN: 209_sex_subtype_figures.R failed (rc=$?)"

echo ""
echo "--- 3. RNA-seq/216_pharma_figures.R (sex panel only relocated) ---"
Rscript "${BASE}/RNA-seq/216_pharma_figures.R" || \
  echo "WARN: 216_pharma_figures.R failed (rc=$?)"

echo ""
echo "============================================================"
echo "  Done. Verifying outputs ..."
echo "============================================================"
ls -la "${BASE}/figures/supplementary/figS_sex_dimorphism/figS_csd_"*.pdf \
       "${BASE}/figures/supplementary/figS_sex_dimorphism/figS_sex_csd_composite.pdf" \
       "${BASE}/figures/supplementary/figS_sex_dimorphism/figS_sex_subtype_"*.pdf \
       "${BASE}/figures/supplementary/figS_sex_dimorphism/figS_sex_pharma_reversal.pdf" 2>&1 | tail -25
