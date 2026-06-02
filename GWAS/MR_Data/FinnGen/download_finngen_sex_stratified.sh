#!/bin/bash -l
# download_finngen_sex_stratified.sh
# ---------------------------------------------------------------------------
# Download FinnGen R12 sex-stratified summary stats for NAFLD/NASH/HCC.
#
# STATUS (2026-05-11): FinnGen R12 PUBLIC release does NOT contain
# sex-stratified summary statistics. See outputs/team_B/B3_blocked.md.
# Sex-stratified GWAS is available only via the FinnGen Sandbox (controlled
# access for approved researchers). This script is staged so it can be
# trivially re-pointed at sandbox-exported files once available.
#
# Once data is available (e.g., manually exported from sandbox to
# /gpfs/commons/.../GWAS/MR_Data/FinnGen/sex_stratified/raw/), set:
#   FINNGEN_SEX_SRC=/path/to/sandbox/exports
# and re-invoke.
#
# Endpoints: NAFLD, CHIRHEP_NAS (= NASH), C3_HEPATOCELLU_CARC_EXALLC (= HCC)
# Per sex: female, male  (6 files total)
# Target dir: GWAS/MR_Data/FinnGen/sex_stratified/
# ---------------------------------------------------------------------------

#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --job-name=fg_sex_dl
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/outputs/team_B/B3_logs/download_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/outputs/team_B/B3_logs/download_%j.err

set -eo pipefail

BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
TARGET_DIR="${BASE_DIR}/GWAS/MR_Data/FinnGen/sex_stratified"
mkdir -p "${TARGET_DIR}"

# Endpoints (FinnGen phenocode -> short label used in registry)
PHENOS_FG=(NAFLD CHIRHEP_NAS C3_HEPATOCELLU_CARC_EXALLC)
PHENOS_LABEL=(NAFLD NASH HCC)
SEXES=(female male)

# Source priority:
#   1. FINNGEN_SEX_SRC env var (manually exported sandbox files)
#   2. Public GCS bucket (will 404 until FinnGen publishes sex-strat)
SRC_DIR="${FINNGEN_SEX_SRC:-}"

echo "=== FinnGen R12 sex-stratified download ==="
echo "Target: ${TARGET_DIR}"
echo "Source override: ${SRC_DIR:-<none, will try public GCS>}"
echo

# If sandbox files were manually exported, just copy
if [[ -n "${SRC_DIR}" && -d "${SRC_DIR}" ]]; then
  echo "Copying from ${SRC_DIR} ..."
  for i in "${!PHENOS_FG[@]}"; do
    pheno_fg="${PHENOS_FG[$i]}"
    pheno_label="${PHENOS_LABEL[$i]}"
    for sex in "${SEXES[@]}"; do
      sex_short=${sex:0:1}      # f / m
      SEX_UPPER=${sex_short^^}  # F / M
      src="${SRC_DIR}/finngen_R12_${pheno_fg}_${sex}.gz"
      dst="${TARGET_DIR}/finngen_R12_${pheno_label}_${SEX_UPPER}.gz"
      if [[ -f "${src}" ]]; then
        cp -v "${src}" "${dst}"
      else
        echo "  MISSING in source: ${src}"
      fi
    done
  done
  exit 0
fi

# Otherwise probe public GCS (will fail today; staged for the day FinnGen publishes)
PUBLIC_BASE="https://storage.googleapis.com/finngen-public-data-r12/summary_stats"
candidates=(
  "${PUBLIC_BASE}/release_female"
  "${PUBLIC_BASE}/release_male"
  "${PUBLIC_BASE}/female"
  "${PUBLIC_BASE}/male"
  "${PUBLIC_BASE}_female"
  "${PUBLIC_BASE}_male"
)

found_any=0
for prefix in "${candidates[@]}"; do
  url="${prefix}/finngen_R12_NAFLD.gz"
  echo -n "Probing ${url} ... "
  code=$(curl -s -o /dev/null -w "%{http_code}" -L "${url}")
  echo "${code}"
  if [[ "${code}" == "200" ]]; then
    found_any=1
    break
  fi
done

if [[ "${found_any}" -eq 0 ]]; then
  echo
  echo "ERROR: FinnGen R12 sex-stratified summary stats not available on public GCS." >&2
  echo "       See outputs/team_B/B3_blocked.md for the documented blocker." >&2
  echo "       To unblock: obtain FinnGen Sandbox access, export the 6 files manually," >&2
  echo "       and re-run with FINNGEN_SEX_SRC=<dir>." >&2
  exit 2
fi

# (Will land here only if FinnGen ever publishes sex-stratified data.)
for i in "${!PHENOS_FG[@]}"; do
  pheno_fg="${PHENOS_FG[$i]}"
  pheno_label="${PHENOS_LABEL[$i]}"
  for sex in "${SEXES[@]}"; do
    sex_short=${sex:0:1}
    SEX_UPPER=${sex_short^^}
    url="${PUBLIC_BASE}/release_${sex}/finngen_R12_${pheno_fg}.gz"
    dst="${TARGET_DIR}/finngen_R12_${pheno_label}_${SEX_UPPER}.gz"
    echo "Downloading ${url}"
    curl -fL --retry 5 --retry-delay 30 -o "${dst}" "${url}"
  done
done

echo "Done."
