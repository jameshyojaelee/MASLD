#!/usr/bin/env bash
#SBATCH --job-name=hac_fig3_validate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/figure-review-validation-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum"

if [[ -z "${HAC_FIGURE_REVIEW_ROOT:-}" ]]; then
  echo "HAC_FIGURE_REVIEW_ROOT is unset" >&2
  exit 1
fi

# The build job already wrote validation.tsv/validation_summary.json into the
# review root, and both are opened write-once. Re-validation gets its own dated
# directory so the original report is never overwritten or blocked.
export HAC_VALIDATION_OUT_DIR="${HAC_VALIDATION_OUT_DIR:-${HAC_FIGURE_REVIEW_ROOT}/revalidation_$(date -u +%Y%m%dT%H%M%SZ)}"
python3 "${SCRIPT_DIR}/11_validate_figure3_review.py"
echo "VALIDATION_REPORT_DIR ${HAC_VALIDATION_OUT_DIR}"
