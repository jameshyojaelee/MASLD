#!/bin/bash
# Generic per-script renderer for the C2 re-render sweep.
# Usage: run_one.sh <chunk_file>  (chunk_file = newline list of script basenames, no .R)
# Writes a per-script result line to $RESDIR/results.tsv : basename<TAB>STATUS<TAB>firstError
# STATUS in {PASS, FAIL, NO_OUTPUT_OK}. Sourced-only detection is post-hoc (parent decides).
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -uo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SDIR=$BASE/scripts/figures
export MASLD_PROJECT_ROOT="$BASE"
RESDIR=$SDIR/_c2_rerender
mkdir -p "$RESDIR/logs"
cd "$BASE"

CHUNK="$1"
RES="$RESDIR/results.tsv"

while read -r name; do
  [ -z "$name" ] && continue
  scr="$SDIR/$name.R"
  log="$RESDIR/logs/$name.log"
  if [ ! -f "$scr" ]; then
    printf '%s\tMISSING\tscript file not found\n' "$name" >> "$RES"
    continue
  fi
  echo "[$(date '+%H:%M:%S')] >>> $name"
  # 20-min per-script cap; most are seconds-to-minutes.
  timeout 1200 Rscript "$scr" > "$log" 2>&1
  rc=$?
  if [ $rc -eq 0 ]; then
    printf '%s\tPASS\t\n' "$name" >> "$RES"
    echo "[$(date '+%H:%M:%S')] <<< $name PASS"
  elif [ $rc -eq 124 ]; then
    printf '%s\tFAIL\tTIMEOUT after 1200s\n' "$name" >> "$RES"
    echo "[$(date '+%H:%M:%S')] <<< $name TIMEOUT"
  else
    # first error-ish line
    ferr=$(grep -m1 -iE "error|cannot|could not|no such file|not found|unable|object .* not found|subscript out of bounds" "$log" | head -1 | tr '\t' ' ' | cut -c1-300)
    [ -z "$ferr" ] && ferr=$(tail -1 "$log" | tr '\t' ' ' | cut -c1-300)
    printf '%s\tFAIL\t%s\n' "$name" "$ferr" >> "$RES"
    echo "[$(date '+%H:%M:%S')] <<< $name FAIL rc=$rc :: $ferr"
  fi
done < "$CHUNK"
echo "[$(date)] chunk $CHUNK DONE"
