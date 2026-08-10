#!/usr/bin/env bash
#SBATCH --job-name=python
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
umask 077
: "${BG001_SOURCE_RUN_ROOT:?BG001_SOURCE_RUN_ROOT is required}"
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"

SOURCE_RUN_ROOT=$(realpath "$BG001_SOURCE_RUN_ROOT")
RUN_ROOT=$(realpath "$BG001_RUN_ROOT")
[[ "$SOURCE_RUN_ROOT" != "$RUN_ROOT" ]] || {
  echo "Source and destination run roots must differ" >&2
  exit 64
}
for root in "$SOURCE_RUN_ROOT" "$RUN_ROOT"; do
  [[ -f "$root/.bg001_candidate_root" && ! -L "$root/.bg001_candidate_root" ]] || {
    echo "Missing or symlinked candidate sentinel: $root" >&2
    exit 64
  }
done

SCRIPT_DIR="$RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
IMPORTER="$SCRIPT_DIR/import_recount.py"
[[ -f "$IMPORTER" && ! -L "$IMPORTER" ]] || {
  echo "Missing frozen destination recount importer: $IMPORTER" >&2
  exit 64
}
mkdir -p "$RUN_ROOT/logs"
exec > >(tee "$RUN_ROOT/logs/import_recount.out") \
  2> >(tee "$RUN_ROOT/logs/import_recount.err" >&2)

python3 "$SCRIPT_DIR/verify_run_contract.py" \
  --run-root "$RUN_ROOT" --require-bam-frozen
python3 "$SCRIPT_DIR/verify_analysis_contract.py" \
  --run-root "$RUN_ROOT" --require-baseline-frozen
python3 "$IMPORTER" \
  --source-run "$SOURCE_RUN_ROOT" \
  --destination-run "$RUN_ROOT" \
  --check-only
python3 "$IMPORTER" \
  --source-run "$SOURCE_RUN_ROOT" \
  --destination-run "$RUN_ROOT"
python3 "$IMPORTER" \
  --verify-import \
  --destination-run "$RUN_ROOT"
