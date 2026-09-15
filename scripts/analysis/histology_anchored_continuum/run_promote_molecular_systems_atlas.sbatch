#!/usr/bin/env bash
#SBATCH --job-name=hac_systems_promote
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROMOTION_SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_molecular_systems_atlas.py"
PROMOTION_MANIFEST="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/manifests/molecular_systems_atlas_promotion_20260824.tsv"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 -m py_compile "${PROMOTION_SCRIPT}"
python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
python3 "${PROMOTION_SCRIPT}"

while IFS=$'\t' read -r callout destination source sha256 rest; do
  [[ "${callout}" == "callout" ]] && continue
  panel="${PROJECT_ROOT}/${destination}"
  [[ "$(sha256sum "${panel}" | awk '{print $1}')" == "${sha256}" ]]
  [[ "$(pdfinfo "${panel}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${panel}"
  ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${panel}"
done < "${PROMOTION_MANIFEST}"

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
printf '%s\n' "promoted_and_validated" > \
  "${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824/PROMOTED.ok"
