#!/bin/bash
# PROT5 step 3: per protein file, extract the +/- 500 kb windows and the null SNPs (data only, zcat + awk).
# Usage: 03_extract.sh <run_dir> <out_dir>
set -euo pipefail
RUN=$1; OUT=$2; X=$OUT/extract; mkdir -p $X
tail -n +2 $RUN/downloads_manifest.tsv | awk -F'\t' '$11=="True"{print $1"\t"$4}' | while IFS=$'\t' read acc f; do
  zcat $RUN/$f | awk -F'\t' -v W=$OUT/windows.tsv -v N=$OUT/null_snps.tsv -v OW=$X/${acc}_win.tsv -v ON=$X/${acc}_null.tsv '
    BEGIN { while ((getline l < W) > 0) { split(l, a, "\t"); n[a[1]]++; s[a[1], n[a[1]]] = a[2] + 0; e[a[1], n[a[1]]] = a[3] + 0 }
            while ((getline l < N) > 0) { split(l, a, "\t"); nul[a[1] ":" a[2]] = 1 } }
    NR == 1 { print > OW; print > ON; next }
    { if (($1 ":" $2) in nul) print > ON
      c = $1; p = $2 + 0
      for (i = 1; i <= n[c]; i++) if (p >= s[c, i] && p <= e[c, i]) { print > OW; break } }'
  gzip -f $X/${acc}_win.tsv $X/${acc}_null.tsv
  echo "$acc extracted: $(zcat $X/${acc}_win.tsv.gz | wc -l) window rows, $(zcat $X/${acc}_null.tsv.gz | wc -l) null rows"
done
