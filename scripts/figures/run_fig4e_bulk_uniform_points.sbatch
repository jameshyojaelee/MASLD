#!/usr/bin/env bash
#SBATCH --job-name=fig4e_uniform
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4e_uniform_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4e_uniform_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export FIG4E_UNIFORM_OUT="${PROJECT_ROOT}/figures/candidates/fig4e-uniform-points-2026-08-17-v1"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

micromamba run -n rnaseq Rscript --vanilla \
  "${PROJECT_ROOT}/scripts/figures/render_fig4e_bulk_uniform_points.R"

pdf="${FIG4E_UNIFORM_OUT}/panels/fig4e_bulk_tissue_state_transport.pdf"
test -s "${pdf}"
test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
if grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${pdf}"; then
  echo "Type 3 font detected: ${pdf}" >&2
  exit 1
fi

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
find "${FIG4E_UNIFORM_OUT}" -type f ! -name output_manifest.tsv -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  | sed "s#  ${FIG4E_UNIFORM_OUT}/#  #" > "${FIG4E_UNIFORM_OUT}/output_manifest.tsv"
touch "${FIG4E_UNIFORM_OUT}/VALIDATED"
echo "FIG4E_UNIFORM_POINTS_VALIDATION PASS"
