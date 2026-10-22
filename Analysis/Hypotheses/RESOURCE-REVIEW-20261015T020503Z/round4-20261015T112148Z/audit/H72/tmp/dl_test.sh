set -euo pipefail
OUT=$1; H72=$1; mkdir -p $H72/downloads
wget() { local o=""; while [ $# -gt 0 ]; do case $1 in -O) o=$2; shift 2;; *) shift;; esac; done; head -c 10 /dev/zero | gzip > "$o"; }
B=x
declare -A LEN=([30620_raw.female]=454804739 [30620_raw.male]=23)
: > "${OUT}/download_status.tsv"
for f in 30620_raw; do for s in female male; do
  fn=${f}.gwas.imputed_v3.${s}.varorder.tsv.bgz; dst=${H72}/downloads/${fn}; ok=0
  for try in 1 2; do
    { [ "${try}" = 1 ] && [ -s "${dst}" ]; } || wget -q -O "${dst}" "${B}/${fn}" || true
    if [ "$(stat -c %s "${dst}" 2>/dev/null)" = "${LEN[${f}.${s}]}" ] && gzip -t "${dst}" 2>/dev/null; then ok=1; break; fi
  done
  printf '%s\t%s\t%s\n' "${fn}" "$(stat -c %s "${dst}" 2>/dev/null)" "${ok}" >> "${OUT}/download_status.tsv"
done; done
cat ${OUT}/download_status.tsv
