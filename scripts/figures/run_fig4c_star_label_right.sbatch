#!/bin/bash
#SBATCH --job-name=fig4c_starlabel
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/fig4c-star-label-right-2026-08-15-v1/logs/render.%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/fig4c-star-label-right-2026-08-15-v1/logs/render.%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SOURCE_ROOT="${PROJECT_ROOT}/figures/candidates/fig4-stage-terminology-corrected-2026-08-13-v4"
CANDIDATE="${PROJECT_ROOT}/figures/candidates/fig4c-star-label-right-2026-08-15-v1"
PROGRAM_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

test -d "${CANDIDATE}"
test -s "${PROGRAM_ROOT}/hotspot/READY"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export FIGURE_CANDIDATE_ROOT="${CANDIDATE}"
export PROGRAM_RELEASE_ROOT="${PROGRAM_ROOT}/hotspot"
export SINGLECELL_STAGE_METADATA_FILE="${PROGRAM_ROOT}/metadata/donor_metadata_stage_corrected.tsv"
export CCC_TRAJECTORY_FILE="${SOURCE_ROOT}/source_tables/fig4d_communication_trajectories.tsv"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/singlecell_module_heatmap.R
"${RNASEQ_ENV}/bin/Rscript" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "sessionInfo.txt"))'
