#!/usr/bin/env bash
#SBATCH --job-name=hac_fig3_review
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/figure-review-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

if [[ -z "${HAC_SOURCE_CANDIDATE:-}" ]]; then
  echo "HAC_SOURCE_CANDIDATE is unset" >&2
  exit 1
fi
if [[ -z "${HAC_FIGURE_REVIEW_ROOT:-}" ]]; then
  echo "HAC_FIGURE_REVIEW_ROOT is unset" >&2
  exit 1
fi

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e '
parse(file = file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
  "scripts/analysis/histology_anchored_continuum/10_render_figure3_review_suite.R"),
  keep.source = FALSE)
cat("FIGURE_REVIEW_R_PARSE_OK\n")
'
"${R_BIN}" --vanilla "${SCRIPT_DIR}/10_render_figure3_review_suite.R"
python3 "${SCRIPT_DIR}/11_validate_figure3_review.py"

echo "FIGURE3_REVIEW_COMPLETE ${HAC_FIGURE_REVIEW_ROOT}"
