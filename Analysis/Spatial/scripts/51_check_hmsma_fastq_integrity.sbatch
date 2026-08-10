#!/bin/bash
#SBATCH --job-name=gzip
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gzip_integrity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gzip_integrity_%j.out
set -o pipefail
# Integrity + read-count check for the 5 HMSMA samples whose STARsolo run died with
# FASTQ desynchronisation / barcode-length errors. Truncated gzip or an R1/R2 read-count
# mismatch both produce exactly those STAR messages.
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
IN=$BASE/Analysis/Spatial/spaceranger_input/HRA007511
OUT=$BASE/Analysis/Spatial/results/starsolo_visium/HRA007511/_integrity_report.tsv

printf "sample\tfile\tgzip_ok\tn_lines\tn_reads\tread_len\n" > "$OUT"
declare -A READ_COUNTS
FAIL=0
for s in HRA_19 HRA_24 HRA_27 HRA_28 HRA_50; do
  for f in "$IN/$s"/*.fastq.gz; do
    if [ ! -e "$f" ]; then
      echo "MISSING $s expected FASTQ glob"
      FAIL=$((FAIL+1))
      continue
    fi
    if gzip -t "$f" 2>/dev/null; then ok=OK; else ok=CORRUPT; FAIL=$((FAIL+1)); fi
    nl=$(zcat "$f" 2>/dev/null | wc -l)
    nr=$((nl/4))
    rl=$(zcat "$f" 2>/dev/null | head -2 | tail -1 | tr -d '\n' | wc -c)
    [ $((nl % 4)) -eq 0 ] || { echo "NON_MULTIPLE_OF_4 $s $(basename "$f") lines=$nl"; FAIL=$((FAIL+1)); }
    [ "$nr" -gt 0 ] || { echo "ZERO_READS $s $(basename "$f")"; FAIL=$((FAIL+1)); }
    READ_COUNTS["$s:$(basename "$f")"]="$nr"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$(basename "$f")" "$ok" "$nl" "$nr" "$rl" >> "$OUT"
    echo "$s $(basename "$f") gzip=$ok reads=$nr len=$rl"
  done

  # Compare every lane's barcode/UMI R1 with its cDNA R2.  Equal paired-read
  # counts are mandatory; total sample counts alone can hide a truncated lane.
  for r1 in "$IN/$s"/*_R1_001.fastq.gz; do
    [ -e "$r1" ] || { echo "MISSING_R1 $s"; FAIL=$((FAIL+1)); continue; }
    r2=${r1/_R1_/_R2_}
    if [ ! -e "$r2" ]; then
      echo "MISSING_MATE $s $(basename "$r1")"
      FAIL=$((FAIL+1))
      continue
    fi
    r1_name=$(basename "$r1")
    r2_name=$(basename "$r2")
    r1_reads=${READ_COUNTS["$s:$r1_name"]:-missing}
    r2_reads=${READ_COUNTS["$s:$r2_name"]:-missing}
    if [ "$r1_reads" = "missing" ] || [ "$r2_reads" = "missing" ] || [ "$r1_reads" != "$r2_reads" ]; then
      echo "PAIR_COUNT_MISMATCH $s $r1_name=$r1_reads $r2_name=$r2_reads"
      FAIL=$((FAIL+1))
    else
      echo "PAIR_COUNT_OK $s $r1_name/$r2_name reads=$r1_reads"
    fi
  done
done
echo "REPORT: $OUT"
if [ "$FAIL" -ne 0 ]; then
  echo "INTEGRITY_FAILED failures=$FAIL"
  exit 1
fi

echo "INTEGRITY_PASSED samples=5"
if [ -n "${STAR_JOB_ID:-}" ]; then
  [[ "$STAR_JOB_ID" =~ ^[0-9]+$ ]] || { echo "INVALID_STAR_JOB_ID=$STAR_JOB_ID"; exit 1; }
  state=$(scontrol show job "$STAR_JOB_ID" | awk -F= '/JobState=/{split($2,a," "); print a[1]; exit}')
  [ "$state" = "PENDING" ] || { echo "STAR_JOB_NOT_PENDING id=$STAR_JOB_ID state=$state"; exit 1; }
  scontrol release "$STAR_JOB_ID"
  echo "STAR_JOB_RELEASED id=$STAR_JOB_ID"
fi
