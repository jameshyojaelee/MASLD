#!/usr/bin/env bash
#SBATCH --job-name=hac_systems_polish_promote
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 "${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_refined_molecular_systems_family.py"

MANIFEST="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/manifests/molecular_systems_visual_revision_promotion_20260824.tsv"
while IFS=$'\t' read -r callout destination source old_sha256 new_sha256 rest; do
  [[ "${callout}" == "callout" ]] && continue
  panel="${PROJECT_ROOT}/${destination}"
  [[ "$(sha256sum "${panel}" | cut -d' ' -f1)" == "${new_sha256}" ]]
  [[ "$(pdfinfo "${panel}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${panel}"
  ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${panel}"
done < "${MANIFEST}"

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
python3 "${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/validate_integration.py"
