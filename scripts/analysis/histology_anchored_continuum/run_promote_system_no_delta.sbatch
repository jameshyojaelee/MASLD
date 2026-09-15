#!/usr/bin/env bash
#SBATCH --job-name=hac_system_no_delta_promote
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 "${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_system_no_delta.py"

MANIFEST="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/manifests/molecular_systems_no_delta_promotion_20260825.tsv"
PANEL="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/panels/supplementary/figs4u_molecular_systems_pooled_continuum_windows.pdf"
EXPECTED="$(awk -F'\t' 'NR==2 {print $5}' "${MANIFEST}")"
[[ "$(sha256sum "${PANEL}" | cut -d' ' -f1)" == "${EXPECTED}" ]]
[[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"
pdftotext "${PANEL}" - | grep -q 'maxT FWER P<0.05'
! pdftotext "${PANEL}" - | grep -q 'Delta'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
