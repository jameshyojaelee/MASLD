#!/bin/bash
#SBATCH --job-name=6121
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/logs/6121-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/logs/6121-%j.err
# Pre-fetch of the eight raw-scale Neale files with the prespec's own download block
# (byte check, gzip -t, one retry). The final job re-runs the same block, which skips
# a file already present on try 1 and re-checks its bytes and gzip integrity.
set -euo pipefail
H72=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z
OUT=${H72}/logs/prefetch-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}"
B=https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs
declare -A LEN=([30620_raw.female]=454804739 [30620_raw.male]=454799513 [30650_raw.female]=454497941
  [30650_raw.male]=454646054 [30730_raw.female]=451369692 [30730_raw.male]=452520192
  [30880_raw.female]=453341535 [30880_raw.male]=454151669)
: > "${OUT}/download_status.tsv"
for f in 30620_raw 30650_raw 30730_raw 30880_raw; do for s in female male; do
  fn=${f}.gwas.imputed_v3.${s}.varorder.tsv.bgz; dst=${H72}/downloads/${fn}; ok=0
  for try in 1 2; do
    { [ "${try}" = 1 ] && [ -s "${dst}" ]; } || wget -q -O "${dst}" "${B}/${fn}" || true
    if [ "$(stat -c %s "${dst}" 2>/dev/null)" = "${LEN[${f}.${s}]}" ] && gzip -t "${dst}" 2>/dev/null; then ok=1; break; fi
  done
  printf '%s\t%s\t%s\n' "${fn}" "$(stat -c %s "${dst}" 2>/dev/null)" "${ok}" >> "${OUT}/download_status.tsv"
done; done
cat "${OUT}/download_status.tsv"
