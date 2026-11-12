#!/bin/bash
# Stream one gzipped GWAS Catalog sumstats file, keep the header and rows inside the hg38 windows, and write the
# sha256 of the compressed file from the same read. Never loads the file into memory.
# Usage: extract_windows.sh <sumstats.tsv.gz> <windows.tsv: chr start end> <sha256_out>
set -euo pipefail
f="$1"; win="$2"; shaout="$3"
tee >(sha256sum | cut -d' ' -f1 > "$shaout") < "$f" | zcat | awk -F'\t' -v W="$win" '
BEGIN { while ((getline l < W) > 0) { split(l, a, "\t"); k = ++n[a[1]]; s[a[1], k] = a[2]; e[a[1], k] = a[3] } }
NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "chromosome") c = i; if ($i == "base_pair_location") p = i }
          if (!c || !p) { print "MISSING_COLUMNS" > "/dev/stderr"; exit 2 }
          print; next }
{ ch = $c; sub(/^chr/, "", ch); if (!(ch in n)) next
  for (j = 1; j <= n[ch]; j++) if ($p + 0 >= s[ch, j] && $p + 0 <= e[ch, j]) { print; break } }'
