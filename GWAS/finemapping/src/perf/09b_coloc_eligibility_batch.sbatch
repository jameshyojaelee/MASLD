#!/bin/bash
#SBATCH --job-name=coloc_elig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_eligibility/coloc_elig_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_eligibility/coloc_elig_%j.err

# COLOC shared-posterior eligibility rerun (review item 2), one batch per job.
# NOT an array: submit once per batch index.
#
#   python3 src/perf/09_prepare_coloc_eligibility_batches.py --run-root <RUN_ROOT>
#   sbatch src/perf/09b_coloc_eligibility_batch.sbatch <RUN_ROOT> <1-5>
#
# Reruns 06_susie_coloc.R (full-fit coloc.susie + per-pair shares) on the
# SuSiE-positive (study, chr, gene) rows of plan/batch_<k>.tsv, with the same
# environment as the adopted rerun (06c_coloc_rerun_full.sbatch): rnaseq env,
# rebuilt Gram LD panels for non-EUR, LD-pinned polyfun eQTL SuSiE fits.
# Writes <RUN_ROOT>/coloc/<gwas>/susie_coloc_chr<N>.csv and
# susie_coloc_pairs_chr<N>.tsv. A (study, chr) task belongs to one batch only.
# Tasks whose two outputs already exist are skipped, so a resubmission after a
# time-out resumes; an interrupted task resumes from 06's own checkpoint.
set -eo pipefail
if [[ $# -ne 2 || ! $2 =~ ^[1-5]$ ]]; then
  echo "Usage: sbatch 09b_coloc_eligibility_batch.sbatch <RUN_ROOT> 1|2|3|4|5" >&2
  exit 2
fi
RUN_ROOT=$1
BATCH=$2
PLAN=${RUN_ROOT}/plan/batch_${BATCH}.tsv
if [[ ! -f "${PLAN}" ]]; then echo "missing ${PLAN}" >&2; exit 2; fi

set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
FM="${PROJECT_ROOT}/GWAS/finemapping"
SRC="${FM}/src"
LDREF="${FM}/data/ld_ref"
export MASLD_PROJECT_ROOT

export AFR_LD_DIR="${LDREF}/1kg_afr_gram"
export AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram"
export SAS_LD_DIR="${LDREF}/1kg_sas_gram"
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
export COLOC_OUT_ROOT="${RUN_ROOT}/coloc"
export COLOC_GENE_FILTER="${PLAN}"

export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# The runner must be the one the plan recorded.
PLANNED_SHA=$(awk -F'\t' '$1=="runner" {print $3}' "${RUN_ROOT}/plan/input_manifest.tsv")
HELPER_SHA=$(awk -F'\t' '$1=="eligibility_helper" {print $3}' "${RUN_ROOT}/plan/input_manifest.tsv")
if [[ "$(sha256sum "${SRC}/06_susie_coloc.R" | cut -d' ' -f1)" != "${PLANNED_SHA}" ||
      "$(sha256sum "${SRC}/coloc_signal_eligibility.R" | cut -d' ' -f1)" != "${HELPER_SHA}" ]]; then
  echo "06_susie_coloc.R or coloc_signal_eligibility.R changed since planning" >&2
  exit 1
fi

mkdir -p "${RUN_ROOT}/env"
Rscript -e 'library(coloc); library(susieR); library(data.table); sessionInfo()' \
  > "${RUN_ROOT}/env/batch_${BATCH}_sessionInfo.txt" 2>&1

cd "${SRC}"
tail -n +2 "${PLAN}" | cut -f2,3 | sort -u | while IFS=$'\t' read -r GWAS CHR; do
  if [[ -f "${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_chr${CHR}.csv" &&
        -f "${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_pairs_chr${CHR}.tsv" ]]; then
    echo "[$(date -u +%FT%TZ)] skip ${GWAS} chr${CHR}: outputs exist"
    continue
  fi
  echo "[$(date -u +%FT%TZ)] batch ${BATCH}: ${GWAS} chr${CHR}"
  GWAS_NAME="${GWAS}" CHR_FILTER="${CHR}" Rscript "${SRC}/06_susie_coloc.R"
done
echo "[$(date -u +%FT%TZ)] batch ${BATCH} done"
