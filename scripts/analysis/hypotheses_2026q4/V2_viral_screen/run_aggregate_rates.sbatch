#!/bin/bash
#SBATCH --job-name=V2agg
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/V2agg/%x-%j.out
#SBATCH --error=logs/V2agg/%x-%j.err

# V2: aggregate finished screen results into a NEW analysis directory without
# touching the original library outputs or the pilot tables. The analysis
# directory holds a tree of symlinks to each library's screen_result.tsv so
# aggregate_screen.py runs unchanged and writes only into the new directory.
#
# Usage: sbatch run_aggregate_rates.sbatch <NEW_ANALYSIS_DIR> [<EXTRA_LIBRARIES_ROOT>]
#   EXTRA_LIBRARIES_ROOT: a libraries/ tree (<cohort>/<library>/screen_result.tsv)
#   from re-screened libraries, added to the tree.

set -euo pipefail

ANALYSIS_DIR="${1:?new analysis dir required}"
EXTRA="${2:-}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
V2="$ROOT/Analysis/Hypotheses/V2-20261007T210714Z"
SCRIPTS="$ROOT/scripts/analysis/hypotheses_2026q4/V2_viral_screen"
IDENT="$ROOT/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

test ! -e "$ANALYSIS_DIR" || { echo "refusing to reuse $ANALYSIS_DIR" >&2; exit 2; }
mkdir -p "$ANALYSIS_DIR/libraries"

link_tree() {
  local src_root="$1" p d lib coh
  for p in "$src_root"/*/*/screen_result.tsv; do
    d=$(dirname "$p"); lib=$(basename "$d"); coh=$(basename "$(dirname "$d")")
    mkdir -p "$ANALYSIS_DIR/libraries/$coh"
    if [ -e "$ANALYSIS_DIR/libraries/$coh/$lib" ]; then
      echo "duplicate library in tree: $coh/$lib" >&2; exit 3
    fi
    ln -s "$d" "$ANALYSIS_DIR/libraries/$coh/$lib"
  done
}
link_tree "$V2/libraries"
if [ -n "$EXTRA" ]; then link_tree "$EXTRA"; fi
echo "libraries linked: $(find "$ANALYSIS_DIR/libraries" -mindepth 2 -maxdepth 2 -type l | wc -l)"

"$PY" -I "$SCRIPTS/aggregate_screen.py" --result-dir "$ANALYSIS_DIR"
"$PY" -I "$SCRIPTS/v2_person_rates.py" --analysis-dir "$ANALYSIS_DIR" \
  --manifest "$V2/library_manifest.tsv" \
  --roster-all "$IDENT/all_screened_roster.tsv" \
  --roster-774 "$IDENT/bulk_execution_inputs/primary_bulk_metadata.tsv" \
  --aggregator "$SCRIPTS/aggregate_screen.py"
"$PY" -I -m pip freeze > "$ANALYSIS_DIR/pip_freeze.txt" 2>&1 || true
"$PY" -I --version > "$ANALYSIS_DIR/python_version.txt" 2>&1
sha256sum "$SCRIPTS/aggregate_screen.py" "$SCRIPTS/score_library.py" "$SCRIPTS/v2_person_rates.py" > "$ANALYSIS_DIR/script_sha256.txt"
echo "done $(date -u +%FT%TZ)"
