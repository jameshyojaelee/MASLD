#!/usr/bin/env bash
#SBATCH --job-name=hac_promote_fig4h
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/slurm-%x-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROMOTER="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_fig4h_system_stage_continuum.py"
PANEL="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/panels/fig4h_molecular_systems_stage_vs_continuum.pdf"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 -m py_compile "${PROMOTER}"
python3 "${PROMOTER}"

[[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"
[[ "$(sha256sum "${PANEL}" | awk '{print $1}')" == \
  "c6098c550f496d8afe427ed7eb0a76d73be3f3f0ca268b4da0f683c2df777615" ]]

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
