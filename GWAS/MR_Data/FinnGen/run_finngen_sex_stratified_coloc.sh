#!/bin/bash -l
# run_finngen_sex_stratified_coloc.sh
# ---------------------------------------------------------------------------
# Submit 6 SuSiE-COLOC array jobs for FinnGen R12 sex-stratified
# (3 endpoints x 2 sexes), staggered 60s apart per CLAUDE.md SLURM protocol.
#
# Prereq: FinnGen sex-stratified gz files in
#   GWAS/MR_Data/FinnGen/sex_stratified/finngen_R12_{NAFLD,NASH,HCC}_{F,M}.gz
# AND registry rows added by format_finngen_sex_stratified.R.
#
# This script is a no-op when sex-strat data is unavailable
# (see outputs/team_B/B3_blocked.md).
#
# Usage: bash run_finngen_sex_stratified_coloc.sh
# ---------------------------------------------------------------------------

set -eo pipefail

BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SEX_DIR="${BASE_DIR}/GWAS/MR_Data/FinnGen/sex_stratified"
REGISTRY="${BASE_DIR}/GWAS/finemapping/config/gwas_registry.tsv"
LOG_DIR="${BASE_DIR}/outputs/team_B/B3_logs"
mkdir -p "${LOG_DIR}"

# Sanity: check inputs exist
all_present=1
for label in NAFLD NASH HCC; do
  for sex in F M; do
    if [[ ! -f "${SEX_DIR}/finngen_R12_${label}_${sex}.gz" ]]; then
      echo "MISSING: ${SEX_DIR}/finngen_R12_${label}_${sex}.gz"
      all_present=0
    fi
  done
done

if [[ "${all_present}" -eq 0 ]]; then
  echo "ERROR: FinnGen sex-stratified gz files missing." >&2
  echo "       See outputs/team_B/B3_blocked.md for the documented blocker." >&2
  echo "       This is expected today; script is staged for future use." >&2
  exit 2
fi

# Sanity: check registry rows
for label in NAFLD NASH HCC; do
  for sex in F M; do
    study="FINNGEN_${sex}_${label}"
    if ! grep -q "^${study}\b" "${REGISTRY}"; then
      echo "MISSING registry row: ${study}" >&2
      echo "Run: Rscript ${BASE_DIR}/GWAS/MR_Data/FinnGen/format_finngen_sex_stratified.R" >&2
      exit 2
    fi
  done
done

# Submit 6 jobs with 60s stagger
cd "${BASE_DIR}/GWAS/finemapping"
JOBS=()
i=0
for label in NAFLD NASH HCC; do
  for sex in F M; do
    study="FINNGEN_${sex}_${label}"
    if [[ ${i} -gt 0 ]]; then
      echo "  Sleeping 60s before next submission..."
      sleep 60
    fi
    echo "Submitting ${study} (array 1-22)..."
    JID=$(GWAS_NAME="${study}" N_WORKERS=4 \
          sbatch --parsable --array=1-22 \
                 --output="${LOG_DIR}/coloc_${study}_%A_%a.out" \
                 --error="${LOG_DIR}/coloc_${study}_%A_%a.err" \
                 src/06_susie_coloc_bigmem.sh)
    echo "  -> Job ID: ${JID}"
    JOBS+=("${JID}:${study}")
    i=$((i + 1))
  done
done

echo
echo "Submitted ${#JOBS[@]} jobs:"
for j in "${JOBS[@]}"; do echo "  ${j}"; done

# Chain aggregation (afterok on all 6 array jobs)
DEPLIST=$(IFS=:; echo "${JOBS[*]%%:*}")
AGG_JID=$(sbatch --parsable \
                 --dependency=afterok:"${DEPLIST}" \
                 --partition=cpu --qos=nslab \
                 --cpus-per-task=8 --mem=32G --time=48:00:00 \
                 --job-name=fg_sex_xtab \
                 --output="${LOG_DIR}/xtab_%j.out" \
                 --error="${LOG_DIR}/xtab_%j.err" \
                 --wrap="cd ${BASE_DIR} && \
                         eval \"\$(micromamba shell hook -s bash)\" && \
                         micromamba activate rnaseq && \
                         Rscript ${BASE_DIR}/RNA-seq/scripts/build_sex_ancestry_coloc_xtab.R")
echo
echo "Aggregation job: ${AGG_JID} (depends on all 6 above)"
