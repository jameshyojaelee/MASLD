#!/bin/bash
#SBATCH --job-name=Rscript
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/rerender_zonation_sensitivity_%j.out
#SBATCH --error=scripts/figures/logs/rerender_zonation_sensitivity_%j.err

# Re-render zonation/sensitivity figures after Conserved-Core (999->1355) +
# zonation PP-in-CC OR (6.20->13.15) recompute (2026-06-10).
# Only re-render the figures whose CONTENT depends on the changed quantities:
#   - figS_lfc_sensitivity.R     (main atlas is_conserved / triple-convergence -> IMPACTED)
#   - figS_spatial_validation.R  (user-listed; spatial replication panels)
# fig4f_{fads2,cyp3a4}_zonation.R are per-gene panels reading the un-recomputed
# spatial-integration atlas; NOT-IMPACTED -> intentionally skipped (see report).

# NOTE: no `set -u` — micromamba's binutils activate.d script references
# unbound vars (ADDR2LINE) and would abort activation under nounset.
set -o pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

cd "${SCRIPT_DIR}" || { echo "FATAL: cannot cd to ${SCRIPT_DIR}"; exit 1; }

echo "[$(date '+%F %T')] host=$(hostname) cwd=$(pwd)"

echo "=========================================="
echo "[$(date '+%T')] figS_lfc_sensitivity.R"
echo "=========================================="
Rscript figS_lfc_sensitivity.R || echo "WARN: figS_lfc_sensitivity.R returned nonzero"

echo "=========================================="
echo "[$(date '+%T')] figS_spatial_validation.R"
echo "=========================================="
Rscript figS_spatial_validation.R || echo "WARN: figS_spatial_validation.R returned nonzero"

echo "[$(date '+%F %T')] DONE"
