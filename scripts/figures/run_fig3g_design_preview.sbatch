#!/bin/bash
#SBATCH --job-name=fig3g_preview
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3g_preview_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3g_preview_%j.err

set -euo pipefail
: "${FIGURE_CANDIDATE_ROOT:?Set fresh preview directory}"
test -d "${FIGURE_CANDIDATE_ROOT}"
test ! -e "${FIGURE_CANDIDATE_ROOT}/DESIGN_PREVIEW_ONLY.txt"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export FIG3G_DESIGN_PREVIEW=true
export PATH="${RNASEQ_ENV}/bin:${PATH}"
"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/figures/fig3g_genetic_state_complementarity.R"
"${RNASEQ_ENV}/bin/Rscript" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "sessionInfo.txt"))'
