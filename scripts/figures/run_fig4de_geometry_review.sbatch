#!/usr/bin/env bash
#SBATCH --job-name=fig4de_geometry
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4de_geometry_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4de_geometry_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export FIG4DE_GEOMETRY_OUT="${PROJECT_ROOT}/figures/candidates/fig4de-geometry-review-2026-08-17-v5"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

micromamba run -n rnaseq Rscript --vanilla \
  "${PROJECT_ROOT}/scripts/figures/render_fig4de_geometry_review.R"

pdfs=(
  "${FIG4DE_GEOMETRY_OUT}/panels/fig4d_lineage_enrichment_forest.pdf"
  "${FIG4DE_GEOMETRY_OUT}/panels/fig4e_independent_stage_contrasts.pdf"
)
for pdf in "${pdfs[@]}"; do
  test -s "${pdf}"
  test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
  if grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${pdf}"; then
    echo "Type 3 font detected: ${pdf}" >&2
    exit 1
  fi
done

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"

find "${FIG4DE_GEOMETRY_OUT}" -type f ! -name output_manifest.tsv -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  | sed "s#  ${FIG4DE_GEOMETRY_OUT}/#  #" > "${FIG4DE_GEOMETRY_OUT}/output_manifest.tsv"
touch "${FIG4DE_GEOMETRY_OUT}/VALIDATED"

echo "FIG4DE_GEOMETRY_REVIEW_VALIDATION PASS"
