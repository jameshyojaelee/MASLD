#!/usr/bin/env bash
#SBATCH --job-name=hac_systems_pool_promote
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROMOTION_SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_pooled_molecular_systems.py"
PROMOTION_MANIFEST="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/manifests/molecular_systems_pooled_promotion_20260824.tsv"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 -m py_compile "${PROMOTION_SCRIPT}"
python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
python3 "${PROMOTION_SCRIPT}"

PANEL="${PROJECT_ROOT}/$(awk -F $'\t' 'NR==2 {print $2}' "${PROMOTION_MANIFEST}")"
EXPECTED="$(awk -F $'\t' 'NR==2 {print $4}' "${PROMOTION_MANIFEST}")"
[[ "$(sha256sum "${PANEL}" | awk '{print $1}')" == "${EXPECTED}" ]]
[[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
python3 "${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/validate_integration.py"
printf '%s\n' "promoted_and_validated" > \
  "${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_pooled_20260824/PROMOTED.ok"
