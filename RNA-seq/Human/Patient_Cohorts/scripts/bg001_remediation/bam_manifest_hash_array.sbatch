#!/usr/bin/env bash
#SBATCH --job-name=samtools
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --array=0-7%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
TASK=${SLURM_ARRAY_TASK_ID:?}
mkdir -p "$BG001_RUN_ROOT/logs"
exec > >(tee "$BG001_RUN_ROOT/logs/bam_hash_${TASK}.out") 2> >(tee "$BG001_RUN_ROOT/logs/bam_hash_${TASK}.err" >&2)

eval "$(micromamba shell hook --shell=bash)"
set +u
micromamba activate rnaseq
set -u
exec bash "$SCRIPT_DIR/hash_bam_manifest_task.sh" --task "$TASK"
