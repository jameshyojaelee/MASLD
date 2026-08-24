#!/usr/bin/env bash
#SBATCH --job-name=fig4def_story
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4def_story_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4def_story_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT="${FIG4DEF_CANDIDATE_ROOT:-${BASE}/figures/candidates/fig4def-multicellular-story-2026-08-17-v2}"
SCRIPT="${BASE}/scripts/figures/render_fig4def_multicellular_story.R"

test ! -e "${OUT}"
export MASLD_PROJECT_ROOT="${BASE}"
export FIG4DEF_CANDIDATE_ROOT="${OUT}"

micromamba run -n rnaseq Rscript "${SCRIPT}"

pdfs=(
  "${OUT}/panels/fig4d_ambient_effect_transport.pdf"
  "${OUT}/panels/fig4e_same_atlas_lineage_specificity.pdf"
  "${OUT}/panels/fig4f_bulk_tissue_state_transport.pdf"
)
for pdf in "${pdfs[@]}"; do
  test -s "${pdf}"
  test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
  if grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${pdf}"; then
    echo "Type 3 font detected: ${pdf}" >&2
    exit 1
  fi
  if gs -dNOSAFER -dNODISPLAY -dPDFDEBUG -f "${pdf}" -c quit 2>&1 \
      | grep -Eq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3'; then
    echo "Type 3 font detected by Ghostscript: ${pdf}" >&2
    exit 1
  fi
done

python3 "${BASE}/scripts/manuscript/validate_resource_scope.py"

find "${OUT}" -type f ! -name output_manifest.tsv -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  | sed "s#  ${OUT}/#  #" > "${OUT}/output_manifest.tsv"
touch "${OUT}/VALIDATED"

echo "FIG4DEF_CANDIDATE_VALIDATION PASS"
