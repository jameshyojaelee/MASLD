#!/bin/bash
set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/Sveinbjornsson2022
TOKEN=06828a13-b8b6-46b7-bfe1-62d1b4eafb8e
BASE="https://download.decode.is/s3/download?token=${TOKEN}&file="
FILES=(NAFL_deCODE_sumstat.txt NAFL_deCODE_sumstat.tmd5sum \
       NAFL_INTERMOUNTAIN_sumstat.txt NAFL_INTERMOUNTAIN_sumstat.tmd5sum \
       NAFL_UKBB_sumstat.txt NAFL_UKBB_sumstat.tmd5sum)
for f in "${FILES[@]}"; do
  echo "=== Downloading $f ==="
  curl -fSL --retry 3 --connect-timeout 30 -o "$f" "${BASE}${f}" || echo "FAILED: $f (token expired?)"
  ls -la "$f" 2>/dev/null || true
done
echo "=== md5 verification ==="
for c in NAFL_deCODE_sumstat NAFL_INTERMOUNTAIN_sumstat NAFL_UKBB_sumstat; do
  if [ -f "${c}.txt" ] && [ -f "${c}.tmd5sum" ]; then
    exp=$(awk 'NR==1{print $1}' "${c}.tmd5sum"); got=$(md5sum "${c}.txt" | awk '{print $1}')
    [ "$exp" = "$got" ] && echo "OK  md5 $c" || echo "MISMATCH $c: exp=$exp got=$got"
  fi
done
echo "=== headers + first data row (column format) ==="
for c in NAFL_deCODE_sumstat NAFL_INTERMOUNTAIN_sumstat NAFL_UKBB_sumstat; do
  [ -f "${c}.txt" ] && { echo "--- $c ---"; head -2 "${c}.txt"; echo "  rows: $(wc -l < ${c}.txt)"; }
done
echo "Done at $(date)"
