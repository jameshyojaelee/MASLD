#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=B6_chrx_master
#SBATCH --output=GWAS/finemapping/logs/B6_chrx_master_%j.out
#SBATCH --error=GWAS/finemapping/logs/B6_chrx_master_%j.err

# ---------------------------------------------------------------------------
# B6 — chrX × MASLD COLOC, master driver.
#
# This single sbatch is the "all-in-one" execution: 58a (extract GTEx chrX
# eQTL + liftover) -> 58 (chrX ABF COLOC for each chrX-bearing GWAS) -> 219
# (cross-tab with A9 chrX classification).
#
# ABF is canonical per C2 (outputs/team_C/C2_critique_sex_inference.md Issue 4).
# SuSiE is NOT run on chrX (male-haploid LD assumption violation).
#
# eQTL substrate: GTEx v8 Liver (N=208). Broadaway has no chrX coverage.
#
# Pooled-sex GWAS (this submission): 2023_36280732_NAFLD_UKBB_EUR is the
# only canonical GWAS with chrX rows (verified empirically — UKBB liver
# enzymes, BBJ, FinnGen, Pan-UKBB non-EUR all dropped chrX during reformat).
# Sex-stratified Pan-UKBB chrX will be added once B1's Neale Round 2
# download lands and we re-extract chrX from the per-sex .varorder.tsv.bgz
# files (separate follow-up).
# ---------------------------------------------------------------------------

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p GWAS/finemapping/logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=========================================================="
echo "B6 chrX × MASLD COLOC master driver"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "=========================================================="

# ---- Step 1: extract GTEx Liver chrX eQTL + liftover -----------------------
EQTL_OUT="data/gtex_v8_liver_eqtl/chrX_marginal_summary_results.tsv"
if [ ! -f "${EQTL_OUT}" ] || [ ! -s "${EQTL_OUT}" ]; then
  echo ""
  echo "[Step 1/3] Extract GTEx v8 Liver chrX eQTL + liftover hg38 -> hg19"
  Rscript GWAS/finemapping/src/58a_extract_gtex_chrx_eqtl.R
else
  echo ""
  echo "[Step 1/3] GTEx chrX eQTL already extracted: ${EQTL_OUT}"
  ls -lh "${EQTL_OUT}"
fi

# ---- Step 2: chrX ABF COLOC for each chrX-bearing GWAS ---------------------
echo ""
echo "[Step 2/3] chrX ABF COLOC (canonical method per C2)"

# Discover GWAS that have chrX rows
CHRX_GWAS=()
REGISTRY_PATH="GWAS/finemapping/config/gwas_registry.tsv"
while IFS=$'\t' read -r study_name sumstats_path leadsnps_path ancestry trait_type n_tot n_cases ld_panel window_mb; do
  if [ "${study_name}" = "study_name" ] || [ -z "${study_name}" ]; then continue; fi
  full="GWAS/finemapping/${sumstats_path}"
  if [ -f "${full}" ]; then
    # Count chrX rows (X or 23)
    cx=$(awk -F'\t' 'NR>1 && ($1=="X" || $1=="23") {c++} END {print c+0}' "${full}")
    if [ "${cx}" -gt 0 ]; then
      echo "   chrX rows in ${study_name}: ${cx}"
      CHRX_GWAS+=("${study_name}")
    fi
  fi
done < "${REGISTRY_PATH}"

if [ ${#CHRX_GWAS[@]} -eq 0 ]; then
  echo "ERROR: No registry GWAS contain chrX rows. Aborting B6."
  exit 1
fi

echo "   GWAS with chrX coverage: ${#CHRX_GWAS[@]}"

# Run COLOC for each
for gw in "${CHRX_GWAS[@]}"; do
  echo ""
  echo "   ---- ${gw} ----"
  GWAS_NAME="${gw}" Rscript GWAS/finemapping/src/58_chrx_coloc_abf.R "${gw}" \
    2>&1 | tee GWAS/finemapping/logs/B6_${gw}_chrx_abf.log
done

# ---- Step 3: cross-tab with A9 ---------------------------------------------
echo ""
echo "[Step 3/3] Cross-tab chrX × sex × COLOC"
Rscript RNA-seq/219_chrx_sex_coloc_xtab.R 2>&1 | tee GWAS/finemapping/logs/B6_xtab.log

echo ""
echo "=========================================================="
echo "B6 chrX × MASLD COLOC complete"
echo "End: $(date)"
echo "=========================================================="
