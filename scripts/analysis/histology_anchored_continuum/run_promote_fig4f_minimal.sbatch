#!/usr/bin/env bash
#SBATCH --job-name=hac_promote_fig4f_minimal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/slurm-%x-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROMOTER="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_fig4f_minimal.py"
PANEL="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/panels/fig4f_continuum_program_trajectories.pdf"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 -m py_compile "${PROMOTER}"
python3 "${PROMOTER}"

[[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"
! pdftotext "${PANEL}" - | grep -Eq 'Source-overlap|Cohort PC1|Fixed projection|Participants'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
