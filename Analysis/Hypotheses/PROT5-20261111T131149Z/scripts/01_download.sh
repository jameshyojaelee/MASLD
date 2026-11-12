#!/bin/bash
# PROT5 step 1: sequential download of the panel pQTL files from the GWAS Catalog FTP, 25 GB cap.
# Usage: 01_download.sh <run_dir>
set -euo pipefail
RUN=$1; DL=$RUN/downloads; MAN=$RUN/downloads_manifest.tsv
CAP=25000000000
BASE=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
mkdir -p $DL
[[ -e $MAN ]] || printf "accession\tgene\tpanel\tfile\turl\tutc_time\tbytes\tsha256\tprovider_md5\tlocal_md5\tmd5_ok\n" > $MAN
total=$(awk -F'\t' 'NR>1{s+=$7} END{print s+0}' $MAN)
tail -n +2 $RUN/inputs/panels.tsv | awk -F'\t' '$1!="retry_not_in_panel" && $4=="True"{print $3"\t"$2"\t"$1}' | while IFS=$'\t' read acc gene panel; do
  if grep -q -P "^$acc\t.*\tTrue$" $MAN; then echo "$acc already verified"; continue; fi
  n=${acc#GCST}; lo=$(( (10#$n - 1) / 1000 * 1000 + 1 )); hi=$(( lo + 999 ))
  dir=$(printf "GCST%08d-GCST%08d/%s" $lo $hi $acc)
  url=$BASE/$dir/$acc.tsv.gz
  size=$(curl -sIL -m 60 "$url" | awk 'tolower($1)=="content-length:"{v=$2} END{gsub("\r","",v); print v+0}')
  if (( size == 0 )); then echo "no Content-Length for $acc; skipping"; continue; fi
  if (( total + size > CAP )); then echo "cap reached before $acc (total $total + $size > $CAP); stopping"; break; fi
  mkdir -p $DL/$acc
  t=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  wget -q -O $DL/$acc/$acc.tsv.gz "$url"
  wget -q -O $DL/$acc/$acc.tsv.gz-meta.yaml "$BASE/$dir/$acc.tsv.gz-meta.yaml" || true
  wget -q -O $DL/$acc/md5sum.txt "$BASE/$dir/md5sum.txt" || true
  b=$(stat -c %s $DL/$acc/$acc.tsv.gz); total=$(( total + b ))
  sha=$(sha256sum $DL/$acc/$acc.tsv.gz | cut -d' ' -f1); md5=$(md5sum $DL/$acc/$acc.tsv.gz | cut -d' ' -f1)
  pmd5=$(awk '/data_file_md5sum/{print $2}' $DL/$acc/$acc.tsv.gz-meta.yaml 2>/dev/null || true)
  [[ -n "$pmd5" ]] || pmd5=$(awk -v f=$acc.tsv.gz '$2==f{print $1}' $DL/$acc/md5sum.txt 2>/dev/null || true)
  ok=False; [[ -n "$pmd5" && "$pmd5" == "$md5" ]] && ok=True
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $acc $gene $panel downloads/$acc/$acc.tsv.gz "$url" $t $b $sha "$pmd5" $md5 $ok >> $MAN
  echo "$acc $gene $b bytes md5_ok=$ok cumulative $total"
done
