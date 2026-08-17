#!/bin/bash
#SBATCH --job-name=pi_fig_finalize
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/pi_fig_finalize_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/pi_fig_finalize_%j.err

set -euo pipefail
umask 027

: "${FIGURE_CANDIDATE_ROOT:?Set the fresh candidate path}"
test -d "${FIGURE_CANDIDATE_ROOT}/analysis/stage_extensions"
test ! -e "${FIGURE_CANDIDATE_ROOT}/output_manifest.tsv"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export PATH="${RNASEQ_ENV}/bin:/gpfs/commons/home/jameslee/.local/bin:${PATH}"

"${RNASEQ_ENV}/bin/python3" "${PROJECT_ROOT}/scripts/figures/finalize_pi_figure_candidate.py"
"${RNASEQ_ENV}/bin/python3" "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
