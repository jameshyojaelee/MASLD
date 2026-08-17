#!/bin/bash
#SBATCH --job-name=fig3d_color
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_color_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_color_%j.err

set -euo pipefail
umask 027

: "${FIG3D_COLOR_CANDIDATE_ROOT:?Set a fresh candidate output directory}"
test -d "${FIG3D_COLOR_CANDIDATE_ROOT}"
test ! -e "${FIG3D_COLOR_CANDIDATE_ROOT}/figure3/panels/fig3d_pca_fibrosis_gradient.pdf"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
unset FIGURE_CANDIDATE_ROOT

"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/figures/figS_pca_definitive.R"
"${RNASEQ_ENV}/bin/Rscript" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIG3D_COLOR_CANDIDATE_ROOT"), "sessionInfo.txt"))'

pdf="${FIG3D_COLOR_CANDIDATE_ROOT}/figure3/panels/fig3d_pca_fibrosis_gradient.pdf"
test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
if strings "${pdf}" | grep -q '/Subtype /Type3'; then
  echo "Type 3 font detected" >&2
  exit 1
fi
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
mkdir -p "${FIG3D_COLOR_CANDIDATE_ROOT}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r200 \
  -sOutputFile="${FIG3D_COLOR_CANDIDATE_ROOT}/proofs/fig3d_pca_fibrosis_gradient.png" \
  "${pdf}"
python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
