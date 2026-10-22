#!/bin/bash
#SBATCH --job-name=V2arr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --array=1-1469%6
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=logs/V2arr-%A_%a.out
#SBATCH --error=logs/V2arr-%A_%a.err

# V2 target screen: full screen over every library in the manifest.
#
# NOT SUBMITTED. This file exists so the exact header can be reviewed. It
# requires explicit approval before it is run.
#
# Usage: sbatch run_screen_array.sbatch <RESULT_DIR>
# <RESULT_DIR>/library_manifest.tsv supplies one library per array index,
# ordered largest FASTQ first.

set -euo pipefail

RESULT_DIR="${1:?result dir required}"
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V2_viral_screen
MANIFEST="$RESULT_DIR/library_manifest.tsv"

read -r IDX COHORT LIBRARY LAYOUT POOLED BYTES FQ1 FQ2 < <(
  awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==i {print; exit}' "$MANIFEST"
)
test -n "${LIBRARY:-}" || { echo "no manifest row for index $SLURM_ARRAY_TASK_ID" >&2; exit 3; }

echo "node=$(hostname) array=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "library=$LIBRARY cohort=$COHORT layout=$LAYOUT bytes_gz=$BYTES"
echo "started=$(date -u +%FT%TZ)"

bash "$SCRIPTS/screen_library.sh" "$LIBRARY" "$COHORT" "$RESULT_DIR" "$FQ1" "${FQ2:-}"

echo "finished=$(date -u +%FT%TZ)"
