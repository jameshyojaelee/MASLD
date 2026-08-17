#!/bin/bash
#SBATCH --job-name=hotspot_stage_fix
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/hotspot_stage_fix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/hotspot_stage_fix_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RELEASE_ID="program-context-v2-stage-corrected-candidate-2026-08-13-v3"
RELEASE_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/${RELEASE_ID}"
RNASEQ_RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

test -d "${RELEASE_ROOT}/hotspot/logs"
test ! -e "${RELEASE_ROOT}/hotspot/gate_status.tsv"
test ! -e "${RELEASE_ROOT}/metadata/donor_metadata_stage_corrected.tsv"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HOTSPOT_V2_RELEASE_ID="${RELEASE_ID}"
export HOTSPOT_V2_METADATA_FILE="${RELEASE_ROOT}/metadata/donor_metadata_stage_corrected.tsv"
export HOTSPOT_POOLED_SCORE_FILE="${RELEASE_ROOT}/hotspot/work/donor_scores_all_weighted.tsv"
export HOTSPOT_V2_SBATCH_NAME="run_hotspot_v2_stage_corrected.sbatch"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_RSCRIPT}" --vanilla \
  Analysis/SingleCell/scripts/hotspot_modules/515_correct_stage_metadata.R
micromamba run -n rnaseq python \
  Analysis/SingleCell/scripts/hotspot_modules/511_pooled_cell_donor_scores.py
"${RNASEQ_RSCRIPT}" --vanilla \
  Analysis/SingleCell/scripts/hotspot_modules/512_hotspot_v2_donor_refit.R
"${RNASEQ_RSCRIPT}" --vanilla \
  Analysis/SingleCell/scripts/hotspot_modules/513_freeze_hotspot_v2_registry.R
"${RNASEQ_RSCRIPT}" --vanilla \
  Analysis/SingleCell/scripts/hotspot_modules/514_validate_hotspot_v2_registry.R

test -s "${RELEASE_ROOT}/hotspot/validation_status.tsv"
test -s "${RELEASE_ROOT}/hotspot/READY"
python3 scripts/manuscript/validate_resource_scope.py
