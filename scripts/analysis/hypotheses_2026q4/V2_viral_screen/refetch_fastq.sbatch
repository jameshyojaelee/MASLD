#!/bin/bash
#SBATCH --job-name=V2refetch
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/V2refetch/%x-%j.out
#SBATCH --error=logs/V2refetch/%x-%j.err

# V2: re-fetch FASTQ for libraries whose local copy is truncated, from ENA only.
# Downloaded files are untrusted data: written to a new empty directory, never
# executed, checked against the ENA-published size and md5, gzip integrity and
# read count before any screening.
#
# Usage: sbatch refetch_fastq.sbatch <NEW_EMPTY_DEST_DIR> <SRR> [<SRR> ...]

set -euo pipefail

DEST="${1:?destination directory required}"
shift
test "$#" -ge 1 || { echo "at least one SRR required" >&2; exit 2; }
if [ -e "$DEST" ] && [ -n "$(ls -A "$DEST" 2>/dev/null)" ]; then
  echo "destination $DEST exists and is not empty" >&2
  exit 2
fi
mkdir -p "$DEST"

echo "node=$(hostname) job=${SLURM_JOB_ID:-na} started=$(date -u +%FT%TZ)"
VER="$DEST/verification.tsv"
printf 'run\tfile\texpected_bytes\tobserved_bytes\texpected_md5\tobserved_md5\tgzip_ok\treads_expected\treads_observed\n' > "$VER"

for SRR in "$@"; do
  REPORT="$DEST/${SRR}.ena_filereport.tsv"
  curl -sS --fail --retry 5 --max-time 120 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=run_accession,library_layout,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" \
    > "$REPORT"
  IFS=$'\t' read -r _run layout reads ftps md5s sizes < <(awk -F'\t' 'NR==2' "$REPORT")
  test "$layout" = "PAIRED" || { echo "$SRR is not PAIRED in ENA: $layout" >&2; exit 3; }
  IFS=';' read -r -a FTP <<< "$ftps"
  IFS=';' read -r -a MD5 <<< "$md5s"
  IFS=';' read -r -a SZ <<< "$sizes"
  test "${#FTP[@]}" -eq 2 || { echo "$SRR: expected 2 FASTQ files, ENA lists ${#FTP[@]}" >&2; exit 3; }
  for i in 0 1; do
    url="https://${FTP[$i]}"
    out="$DEST/$(basename "${FTP[$i]}")"
    echo "download $url -> $out $(date -u +%FT%TZ)"
    curl -sS --fail -L --retry 5 --retry-delay 30 -C - -o "$out" "$url"
    obs_bytes=$(stat -c %s "$out")
    obs_md5=$(md5sum "$out" | cut -d' ' -f1)
    if gzip -t "$out" 2> "$out.gziptest.err"; then gz=yes; else gz=no; fi
    n=$(( $(zcat "$out" | wc -l) / 4 )) || n=NA
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SRR" "$(basename "$out")" "${SZ[$i]}" "$obs_bytes" "${MD5[$i]}" "$obs_md5" "$gz" "$reads" "$n" >> "$VER"
  done
done
echo "finished=$(date -u +%FT%TZ)"
cat "$VER"
