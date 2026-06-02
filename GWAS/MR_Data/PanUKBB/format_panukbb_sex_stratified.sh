#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=panukbb_sex_fmt
#SBATCH --output=GWAS/MR_Data/PanUKBB/sex_stratified/logs/format_%j.out
#SBATCH --error=GWAS/MR_Data/PanUKBB/sex_stratified/logs/format_%j.err

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# rnaseq env's binutils activate hook references unbound ADDR2LINE; relax `-u` around activate.
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

echo "=== Sex-stratified PanUKBB formatting ==="
echo "Start: $(date)"

Rscript GWAS/MR_Data/PanUKBB/format_panukbb_sex_stratified.R

# ---------------------------------------------------------------------------
# Append registry entries (idempotent — uses awk to skip if already present).
# Sample sizes are pulled from the format step's summary table.
# ---------------------------------------------------------------------------
REG=GWAS/finemapping/config/gwas_registry.tsv
SUMMARY=GWAS/MR_Data/PanUKBB/sex_stratified/format_summary.tsv

if [ -s "${SUMMARY}" ]; then
  echo ""
  echo "--- Updating GWAS registry ---"
  while IFS=$'\t' read -r stratum trait sex n_var n_sig n_samp; do
    [ "${stratum}" = "stratum" ] && continue
    if grep -q "^${stratum}	" "${REG}"; then
      echo "  Already present: ${stratum}"
      continue
    fi
    echo "  Appending: ${stratum} (N=${n_samp})"
    printf '%s\tdata/sumstats/%s_reformatted_hg19.tsv\tdata/lead_snps/%s_leadSNPs.tsv\tEUR\tquantitative\t%s\t0\tukbb_eur\t0.5\n' \
      "${stratum}" "${stratum}" "${stratum}" "${n_samp}" >> "${REG}"
  done < "${SUMMARY}"
  echo "  Registry tail:"
  tail -8 "${REG}"
fi

echo ""
echo "End: $(date)"
