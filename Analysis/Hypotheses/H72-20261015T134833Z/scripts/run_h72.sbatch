#!/bin/bash
#SBATCH --job-name=6120
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/slurm-6120-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/slurm-6120-%j.err
set -euo pipefail
H72=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z
OUT=${H72}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
B=https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs
declare -A LEN=([30620_raw.female]=454804739 [30620_raw.male]=454799513 [30650_raw.female]=454497941
  [30650_raw.male]=454646054 [30730_raw.female]=451369692 [30730_raw.male]=452520192
  [30880_raw.female]=453341535 [30880_raw.male]=454151669)
: > "${OUT}/download_status.tsv"
for f in 30620_raw 30650_raw 30730_raw 30880_raw; do for s in female male; do
  fn=${f}.gwas.imputed_v3.${s}.varorder.tsv.bgz; dst=${H72}/downloads/${fn}; ok=0
  for try in 1 2; do   # S2: one retry
    { [ "${try}" = 1 ] && [ -s "${dst}" ]; } || wget -q -O "${dst}" "${B}/${fn}" || true
    if [ "$(stat -c %s "${dst}" 2>/dev/null)" = "${LEN[${f}.${s}]}" ] && gzip -t "${dst}" 2>/dev/null; then ok=1; break; fi
  done
  printf '%s\t%s\t%s\n' "${fn}" "$(stat -c %s "${dst}" 2>/dev/null)" "${ok}" >> "${OUT}/download_status.tsv"
done; done
# 01_h72_sex_effects.py reads download_status.tsv first: any ok=0 applies S2 (enzyme file) or S6 (urate file)
stat -c '%s %n' ${H72}/downloads/*.bgz > "${OUT}/download_bytes.txt"
sha256sum ${H72}/downloads/*.bgz > "${OUT}/download_sha256.txt"
${ENV}/bin/python ${H72}/scripts/01_h72_sex_effects.py --out "${OUT}" --downloads ${H72}/downloads --seed 20260923
${ENV}/bin/Rscript ${H72}/scripts/02_h72_figure.R "${OUT}" 2> "${OUT}/figure_captions.txt"
${ENV}/bin/python -m pip freeze > "${OUT}/pip_freeze.txt"
sha256sum ${H72}/scripts/* > "${OUT}/script_sha256.tsv"
