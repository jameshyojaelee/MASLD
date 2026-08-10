#!/usr/bin/env bash
#SBATCH --job-name=samtools
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
umask 077
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
mkdir -p "$BG001_RUN_ROOT/logs"
exec > >(tee "$BG001_RUN_ROOT/logs/bam_manifest_guarded_finalize.out") \
  2> >(tee "$BG001_RUN_ROOT/logs/bam_manifest_guarded_finalize.err" >&2)

eval "$(micromamba shell hook --shell=bash)"
set +u
micromamba activate rnaseq
set -u

python3 "$SCRIPT_DIR/bam_manifest_finalization_guard.py" pre --run-root "$BG001_RUN_ROOT"

# Execute the byte-identical legacy finalization producer between the two guard
# phases.  Calling it as a child preserves this wrapper's postflight step while
# its embedded SBATCH directives remain inert.
bash "$SCRIPT_DIR/bam_manifest_finalize.sbatch"

python3 "$SCRIPT_DIR/bam_manifest_finalization_guard.py" post --run-root "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/bam_manifest_finalization_guard.py" verify --run-root "$BG001_RUN_ROOT"
