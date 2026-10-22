#!/bin/bash
#SBATCH --job-name=V2sra
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/V2refetch/%x-%j.out
#SBATCH --error=logs/V2refetch/%x-%j.err

# V2: re-fetch one library whose local FASTQ is truncated, from the NCBI SRA
# open-data copy (prefetch, then fasterq-dump). Used instead of the ENA FASTQ
# route (refetch_fastq.sbatch), which ran at about 0.65 MB/s.
#
# Downloaded files are untrusted data: written to a new empty directory, never
# executed, and checked before any screening: vdb-validate on the .sra, gzip
# integrity of each FASTQ, and read counts of both mates equal to the read count
# ENA publishes for the run. The job exits non-zero if any check fails, so a
# dependent screening job (--dependency=afterok) never starts on a bad file.
# Uncompressed fasterq-dump output is kept in <DEST>/fastq_raw.
#
# Usage: sbatch fetch_sra_fastq.sbatch <NEW_EMPTY_DEST_DIR> <SRR>

set -euo pipefail

DEST="${1:?destination directory required}"
SRR="${2:?SRR required}"
if [ -e "$DEST" ] && [ -n "$(ls -A "$DEST" 2>/dev/null)" ]; then
  echo "destination $DEST exists and is not empty" >&2
  exit 2
fi
mkdir -p "$DEST/sra" "$DEST/tmp" "$DEST/fastq_raw" "$DEST/fastq"

module load SRA-Toolkit/3.2.1
PIGZ=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pigz
T="${SLURM_CPUS_PER_TASK:-8}"

echo "node=$(hostname) job=${SLURM_JOB_ID:-na} srr=$SRR started=$(date -u +%FT%TZ)"
{
  echo "prefetch $(prefetch --version 2>&1 | tail -1)"
  echo "fasterq-dump $(fasterq-dump --version 2>&1 | grep -i 'fasterq-dump' | tail -1)"
  echo "pigz $("$PIGZ" --version 2>&1)"
} > "$DEST/tool_versions.txt"

# ---------------------------------------------------- expected values (ENA)
REPORT="$DEST/${SRR}.ena_filereport.tsv"
curl -sS --fail --retry 5 --max-time 120 \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=run_accession,library_layout,read_count,fastq_bytes&format=tsv" \
  > "$REPORT"
IFS=$'\t' read -r _run layout reads _sizes < <(awk -F'\t' 'NR==2' "$REPORT")
test "$layout" = "PAIRED" || { echo "$SRR is not PAIRED in ENA: $layout" >&2; exit 3; }
test -n "$reads" || { echo "no ENA read count for $SRR" >&2; exit 3; }
SRA_URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR}"
exp_sra_bytes=$(curl -sI --max-time 60 "$SRA_URL" | awk 'tolower($1)=="content-length:"{print $2+0}')
echo "ena_read_count=$reads expected_sra_bytes=${exp_sra_bytes:-NA}"

# ------------------------------------------------------------- prefetch
for i in 1 2 3; do
  prefetch "$SRR" -O "$DEST/sra" --max-size 50G && break
  echo "prefetch retry $i"; sleep 60
done
SRA="$DEST/sra/$SRR/$SRR.sra"
test -s "$SRA" || { echo "prefetch produced no $SRA" >&2; exit 4; }
obs_sra_bytes=$(stat -c %s "$SRA")
echo "observed_sra_bytes=$obs_sra_bytes"

vdb-validate "$SRA" > "$DEST/vdb_validate.txt" 2>&1 || true
cat "$DEST/vdb_validate.txt"
vdb=no
if grep -q "is consistent" "$DEST/vdb_validate.txt" && ! grep -qiE "inconsistent|corrupt" "$DEST/vdb_validate.txt"; then
  vdb=yes
fi
test "$vdb" = "yes" || { echo "vdb-validate did not report a consistent file" >&2; exit 5; }

# ------------------------------------------------------------ fasterq-dump
fasterq-dump "$SRA" --split-3 -e "$T" -m 8G -t "$DEST/tmp" -O "$DEST/fastq_raw"
ls -la "$DEST/fastq_raw"

VER="$DEST/verification.tsv"
printf 'run\tfile\tena_read_count\tobserved_reads\tbytes_gz\tgzip_ok\tmd5_gz\tvdb_validate\texpected_sra_bytes\tobserved_sra_bytes\n' > "$VER"
status=0
for m in 1 2; do
  raw="$DEST/fastq_raw/${SRR}_${m}.fastq"
  test -s "$raw" || { echo "missing $raw" >&2; exit 6; }
  nl=$(wc -l < "$raw")
  test $(( nl % 4 )) -eq 0 || { echo "$raw line count $nl is not a multiple of 4" >&2; status=1; }
  obs=$(( nl / 4 ))
  gz="$DEST/fastq/${SRR}_${m}.fastq.gz"
  "$PIGZ" -p "$T" -c "$raw" > "$gz"
  if gzip -t "$gz" 2> "$gz.gziptest.err"; then gzok=yes; else gzok=no; status=1; fi
  [ "$obs" -eq "$reads" ] || status=1
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SRR" "$(basename "$gz")" "$reads" "$obs" \
    "$(stat -c %s "$gz")" "$gzok" "$(md5sum "$gz" | cut -d' ' -f1)" "$vdb" "${exp_sra_bytes:-NA}" "$obs_sra_bytes" >> "$VER"
done
if [ -e "$DEST/fastq_raw/${SRR}.fastq" ]; then
  echo "unpaired reads file present: $(( $(wc -l < "$DEST/fastq_raw/${SRR}.fastq") / 4 )) reads" | tee "$DEST/unpaired_reads.txt"
fi
cat "$VER"
test "$status" -eq 0 || { echo "verification failed" >&2; exit 7; }
echo "FETCH_OK $SRR" | tee "$DEST/FETCH_OK"
echo "finished=$(date -u +%FT%TZ)"
