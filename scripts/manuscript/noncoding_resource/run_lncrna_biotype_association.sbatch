#!/bin/bash
#SBATCH --job-name=ncbiotype
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncbiotype_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncbiotype_%j.err

set -euo pipefail
umask 027

: "${NONCODING_SOURCE_GATE:?Set NONCODING_SOURCE_GATE}"
: "${NONCODING_IDENTITY_ROOT:?Set NONCODING_IDENTITY_ROOT}"
: "${NONCODING_BIOTYPE_OUTPUT:?Set NONCODING_BIOTYPE_OUTPUT}"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/build_lncrna_biotype_association.R" \
  "--source-gate=${NONCODING_SOURCE_GATE}" \
  "--gene-identity=${NONCODING_IDENTITY_ROOT}/gene_identity.tsv" \
  "--output=${NONCODING_BIOTYPE_OUTPUT}"
