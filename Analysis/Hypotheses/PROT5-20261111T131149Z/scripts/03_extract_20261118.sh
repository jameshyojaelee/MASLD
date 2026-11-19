#!/bin/bash
# PROT5 step 3 (2026-11-18 remediation): parallel per-protein extraction of +/- 500 kb windows and null SNPs.
# Usage: 03_extract_20261118.sh <run_dir> <out_dir> [n_jobs]
set -euo pipefail
RUN=$1; OUT=$2; NJOBS=${3:-4}
X=$OUT/extract; mkdir -p "$X"

extract_one() {
  local acc="$1" f="$2"
  zcat "$RUN/$f" | awk -F'\t' -v W="$OUT/windows.tsv" -v N="$OUT/null_snps.tsv" \
    -v OW="$X/${acc}_win.tsv" -v ON="$X/${acc}_null.tsv" '
    BEGIN {
      while ((getline l < W) > 0) { split(l, a, "\t"); n[a[1]]++; s[a[1], n[a[1]]] = a[2] + 0; e[a[1], n[a[1]]] = a[3] + 0 }
      while ((getline l < N) > 0) { split(l, a, "\t"); nul[a[1] ":" a[2]] = 1 }
    }
    NR == 1 { print > OW; print > ON; next }
    {
      if (($1 ":" $2) in nul) print > ON
      c = $1; p = $2 + 0
      for (i = 1; i <= n[c]; i++) if (p >= s[c, i] && p <= e[c, i]) { print > OW; break }
    }'
  gzip -f "$X/${acc}_win.tsv" "$X/${acc}_null.tsv"
  echo "$acc extracted: $(zcat "$X/${acc}_win.tsv.gz" | wc -l) window rows, $(zcat "$X/${acc}_null.tsv.gz" | wc -l) null rows"
}
export -f extract_one
export RUN OUT X

tail -n +2 "$RUN/downloads_manifest.tsv" | awk -F'\t' '$11=="True"{print $1"\t"$4}' \
  | xargs -P "$NJOBS" -n 2 bash -c 'extract_one "$1" "$2"' _
