#!/bin/bash
#SBATCH --job-name=fig4_ccc_stage_fix
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_ccc_stage_fix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_ccc_stage_fix_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROGRAM_RELEASE="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3"
FIGURE_ROOT="${PROJECT_ROOT}/figures/candidates/fig4-stage-terminology-corrected-2026-08-13-v3"
RNASEQ_RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

test -s "${PROGRAM_RELEASE}/hotspot/READY"
test -s "${PROGRAM_RELEASE}/metadata/donor_metadata_stage_corrected.tsv"
test -s "${PROGRAM_RELEASE}/metadata/donor_metadata_stage_corrected_donor.tsv"
test -d "${FIGURE_ROOT}"
test ! -e "${FIGURE_ROOT}/source_tables/fig4d_communication_audit.tsv"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export SINGLECELL_STAGE_METADATA_FILE="${PROGRAM_RELEASE}/metadata/donor_metadata_stage_corrected.tsv"
export SINGLECELL_STAGE_DONOR_METADATA_FILE="${PROGRAM_RELEASE}/metadata/donor_metadata_stage_corrected_donor.tsv"
export CCC_CANDIDATE_ROOT="${FIGURE_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_RSCRIPT}" --vanilla scripts/figures/ccc_v3_dc_rebuild.R
python3 scripts/manuscript/validate_resource_scope.py
