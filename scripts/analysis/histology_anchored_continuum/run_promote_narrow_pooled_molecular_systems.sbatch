#!/usr/bin/env bash
#SBATCH --job-name=hac_systems_pool_narrow_promote
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 "${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_narrow_pooled_molecular_systems.py"

MANIFEST="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/manifests/molecular_systems_pooled_narrow_promotion_20260824.tsv"
IFS=$'\t' read -r callout destination source old_sha256 new_sha256 rest < <(tail -n 1 "${MANIFEST}")
panel="${PROJECT_ROOT}/${destination}"
[[ "${callout}" == "S4U" ]]
[[ "$(sha256sum "${panel}" | cut -d' ' -f1)" == "${new_sha256}" ]]
[[ "$(pdfinfo "${panel}" | awk '/^Pages:/ {print $2}')" == "1" ]]
[[ "$(pdfinfo "${panel}" | awk '/^Page size:/ {print $3 "x" $5}')" == "162x370" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${panel}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${panel}"
pdftotext "${panel}" - | grep -q 'Delta = mean(80,90) - mean(10,20).'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
python3 "${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/validate_integration.py"
