#!/bin/bash
#SBATCH --job-name=fig4ef_expand
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4ef_expand_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4ef_expand_%j.err

set -euo pipefail
umask 027

: "${FIGURE_CANDIDATE_ROOT:?Set a fresh candidate output directory}"
: "${STAGE_RELEASE_ROOT:?Set the validated five-cohort fragment release}"
test -d "${FIGURE_CANDIDATE_ROOT}"
test -d "${FIGURE_CANDIDATE_ROOT}/analysis/stage_extensions"
test ! -e "${FIGURE_CANDIDATE_ROOT}/figure4/panels/fig4e_bulk_projection_pre_genetics.pdf"
test ! -e "${FIGURE_CANDIDATE_ROOT}/figure4/panels/fig4f_tf_activity.pdf"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/figures/render_pi_singlecell_panels.R"
"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/figures/tf_convergence_vsF0.R"
"${RNASEQ_ENV}/bin/Rscript" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "sessionInfo.txt"))'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
