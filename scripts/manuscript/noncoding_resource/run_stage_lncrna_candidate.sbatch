#!/bin/bash
#SBATCH --job-name=ncstage
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncstage_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncstage_%j.err

set -euo pipefail
umask 027

: "${NONCODING_SOURCE_GATE:?Set NONCODING_SOURCE_GATE to source_gate_status.tsv}"
: "${NONCODING_IDENTITY_ROOT:?Set NONCODING_IDENTITY_ROOT to the immutable identity candidate}"
: "${NONCODING_STAGE_OUTPUT:?Set NONCODING_STAGE_OUTPUT to a new output directory}"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
STAGE_PREFLIGHT="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-STAGE/PREFLIGHT-v1"
PRODUCER="${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/build_stage_lncrna_candidate.R"
VALIDATOR="${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/validate_stage_lncrna.py"

test -f "${NONCODING_SOURCE_GATE}"
test -f "${NONCODING_IDENTITY_ROOT}/gene_identity.tsv"
test -f "${NONCODING_IDENTITY_ROOT}/lncrna_genomic_class.tsv"
test -f "${STAGE_PREFLIGHT}/STAGE_PREFLIGHT_READY.json"
test ! -e "${NONCODING_STAGE_OUTPUT}"
test ! -e "${NONCODING_STAGE_OUTPUT}.review"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PATH="${RNASEQ_ENV}/bin:${PATH}"
export CONDA_PREFIX="${RNASEQ_ENV}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla "${PRODUCER}" \
  "--source-gate=${NONCODING_SOURCE_GATE}" \
  "--gene-identity=${NONCODING_IDENTITY_ROOT}/gene_identity.tsv" \
  "--lncrna-class=${NONCODING_IDENTITY_ROOT}/lncrna_genomic_class.tsv" \
  "--stage-preflight=${STAGE_PREFLIGHT}" \
  "--output=${NONCODING_STAGE_OUTPUT}"

"${RNASEQ_ENV}/bin/python3" "${VALIDATOR}" \
  --candidate "${NONCODING_STAGE_OUTPUT}" \
  --output "${NONCODING_STAGE_OUTPUT}.review"

test -f "${NONCODING_STAGE_OUTPUT}.review/VALIDATED.json"
echo "STAGE_LNCRNA_REVIEW_COMPLETE: ${NONCODING_STAGE_OUTPUT}.review/VALIDATED.json"
