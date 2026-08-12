#!/bin/bash
#SBATCH --job-name=ncfig3
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncfig3_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncfig3_%j.err

set -euo pipefail
umask 027

: "${NONCODING_BIOTYPE_ROOT:?Set NONCODING_BIOTYPE_ROOT to the validated pooled candidate}"
: "${NONCODING_BIOTYPE_REVIEW:?Set NONCODING_BIOTYPE_REVIEW to VALIDATED.json}"
: "${NONCODING_STAGE_ROOT:?Set NONCODING_STAGE_ROOT to the validated stage candidate}"
: "${NONCODING_STAGE_REVIEW:?Set NONCODING_STAGE_REVIEW to VALIDATED.json}"
: "${NONCODING_FIGURE_OUTPUT:?Set NONCODING_FIGURE_OUTPUT to a new output directory}"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
PRODUCER="${PROJECT_ROOT}/scripts/figures/noncoding_resource_fig3_panels.R"

test -f "${NONCODING_BIOTYPE_ROOT}/bulk_deg_counts_by_biotype.tsv"
test -f "${NONCODING_BIOTYPE_REVIEW}"
test -f "${NONCODING_STAGE_ROOT}/stage_counts_by_biotype.tsv"
test -f "${NONCODING_STAGE_REVIEW}"
test ! -e "${NONCODING_FIGURE_OUTPUT}"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PATH="${RNASEQ_ENV}/bin:${PATH}"
export CONDA_PREFIX="${RNASEQ_ENV}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla "${PRODUCER}" \
  "--bulk-counts=${NONCODING_BIOTYPE_ROOT}/bulk_deg_counts_by_biotype.tsv" \
  "--bulk-review=${NONCODING_BIOTYPE_REVIEW}" \
  "--stage-counts=${NONCODING_STAGE_ROOT}/stage_counts_by_biotype.tsv" \
  "--stage-review=${NONCODING_STAGE_REVIEW}" \
  "--output=${NONCODING_FIGURE_OUTPUT}"
