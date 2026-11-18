#!/bin/bash
#SBATCH --job-name=stage_meta_pi_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/stage_meta_pi_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/stage_meta_pi_20261118_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
CAND_META="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/candidate_20261118"
mkdir -p "${CAND_META}" "${PROJECT_ROOT}/logs"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HARMONIZE_OUT_DIR="${CAND_META}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"${R_BIN}" --vanilla "${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/00_harmonize_metadata.R"
"${R_BIN}" --vanilla "${PROJECT_ROOT}/scripts/analysis/bulk_refit_followup/prepare_stage_metadata_20261118.R"

export STAGE_RELEASE_ROOT="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/stage-corrected-20261118/primary_source_controls"
export STAGE_CONTEXT_ROSTER="${CAND_META}/primary_bulk_metadata_20261118.tsv"
export STAGE_FIXED_NMF_SCORES="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/fixed-nmf-corrected-scores-20261009T141312Z/fixed_nmf_scores_full_stage_k6.tsv"
export FIGURE_CANDIDATE_ROOT="${PROJECT_ROOT}/figures/candidates/pi-stage-20261118"
rm -rf "${FIGURE_CANDIDATE_ROOT}"

"${R_BIN}" --vanilla "${PROJECT_ROOT}/scripts/figures/build_pi_stage_extensions.R"
echo "COMPLETE stage_meta_pi_20261118"
