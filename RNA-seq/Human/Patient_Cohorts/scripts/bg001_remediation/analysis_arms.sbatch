#!/usr/bin/env bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
mkdir -p "$BG001_RUN_ROOT/logs"
exec > >(tee "$BG001_RUN_ROOT/logs/analysis_arms.out") 2> >(tee "$BG001_RUN_ROOT/logs/analysis_arms.err" >&2)
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT" --require-bam-frozen
exec bash "$SCRIPT_DIR/run_four_arms.sh" "$BG001_RUN_ROOT"
