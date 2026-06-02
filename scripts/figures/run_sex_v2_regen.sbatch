#!/bin/bash
#SBATCH --job-name=sex_v2_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_v2_figs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_v2_figs_%j.err

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="${BASE}/scripts/figures"

mkdir -p "${BASE}/figures/logs"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "============================================================"
echo "  REGENERATE SEX-RELATED FIGURES (v2 interaction-based)"
echo "  Date: $(date)"
echo "  Job ID: ${SLURM_JOB_ID:-local}"
echo "============================================================"

echo ""
echo "--- 1. figS_sex_divergent.R (sex divergent sensitivity) ---"
Rscript "${SCRIPT_DIR}/figS_sex_divergent.R"

echo ""
echo "--- 2. fig2_disease_progression.R (sex scatter + class panels) ---"
Rscript "${SCRIPT_DIR}/fig2_disease_progression.R"

echo ""
echo "--- 3. fig5_convergence.R (evidence matrix with sex column) ---"
Rscript "${SCRIPT_DIR}/fig5_convergence.R"

echo ""
echo "--- 4. figS_evidence_convergence.R (evidence convergence with sex) ---"
Rscript "${SCRIPT_DIR}/figS_evidence_convergence.R"

echo "============================================================"
echo "  SEX V2 FIGURE REGENERATION COMPLETE"
echo "  Finished: $(date)"
echo "============================================================"
echo ""
echo "Outputs:"
ls -lh "${BASE}/figures/supplementary/figS03_deconvolution/"figS*sex* 2>/dev/null || true
ls -lh "${BASE}/figures/main/fig2_f2_switch/"*disease_progression* 2>/dev/null || true
ls -lh "${BASE}/figures/main/fig5_convergence/"*evidence_matrix* 2>/dev/null || true
ls -lh "${BASE}/figures/supplementary/"*convergence* 2>/dev/null || true
