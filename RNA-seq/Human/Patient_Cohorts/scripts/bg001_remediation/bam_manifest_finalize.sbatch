#!/usr/bin/env bash
#SBATCH --job-name=samtools
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
mkdir -p "$BG001_RUN_ROOT/logs"
exec > >(tee "$BG001_RUN_ROOT/logs/bam_manifest_finalize.out") 2> >(tee "$BG001_RUN_ROOT/logs/bam_manifest_finalize.err" >&2)
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/finalize_bam_manifest.py" --run-root "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT" --require-bam-frozen
