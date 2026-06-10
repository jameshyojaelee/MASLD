#!/bin/bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/fig4_revalidation_%j.out
#SBATCH --error=scripts/figures/logs/fig4_revalidation_%j.err

# Re-render validation/proteomics Fig 4 panels after Conserved-Core +
# proteomics-concordance recompute (CC 999->1,355; CC protein rho 0.563->0.519).
# Data already recomputed; this only re-renders. PDF only.

set -eo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

cd "${SCRIPT_DIR}"

# Dependency-aware order:
#   1) module_evidence_scan writes panels/data/module_evidence_matrix.csv
#      (required by fig4_validation_v2.R panel 4a)
#   2) fig4_validation.R writes pxd052937_mrna_protein_concordance.csv
#      (fallback input for fig4_validation_v2.R panel 4d)
#   3) fig4_validation_v2.R (overwrites fig4a-d.pdf in panels/ with v2 design)
#   4) fig4a_proteomics_overview.R (standalone)
SCRIPTS=(
  "fig4_module_evidence_scan.R"
  "fig4_validation.R"
  "fig4_validation_v2.R"
  "fig4a_proteomics_overview.R"
)

for s in "${SCRIPTS[@]}"; do
  echo "=========================================="
  echo "[$(date '+%H:%M:%S')] Running ${s}"
  echo "=========================================="
  Rscript "${SCRIPT_DIR}/${s}" || echo "WARN: ${s} FAILED"
  echo ""
done

echo "[$(date '+%H:%M:%S')] DONE"
ls -la --time-style=full-iso "${BASE}/figures/main/fig4_validation/"*.pdf \
  "${BASE}/figures/main/fig4_validation/panels/"*.pdf 2>/dev/null
