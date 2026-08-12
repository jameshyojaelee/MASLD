#!/bin/bash
#SBATCH --job-name=nclncrna
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclncrna_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclncrna_%j.err

set -euo pipefail
umask 027

: "${NONCODING_SOURCE_GATE:?Set NONCODING_SOURCE_GATE to source_gate_status.tsv}"
: "${NONCODING_IDENTITY_ROOT:?Set NONCODING_IDENTITY_ROOT to the immutable identity candidate}"
: "${NONCODING_BULK_OUTPUT:?Set NONCODING_BULK_OUTPUT to a new output directory}"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
PRODUCER="${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/build_bulk_lncrna_candidate.R"

test -f "${NONCODING_SOURCE_GATE}"
test -f "${NONCODING_IDENTITY_ROOT}/gene_identity.tsv"
test -f "${NONCODING_IDENTITY_ROOT}/lncrna_genomic_class.tsv"
test ! -e "${NONCODING_BULK_OUTPUT}"

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
  "--output=${NONCODING_BULK_OUTPUT}"
