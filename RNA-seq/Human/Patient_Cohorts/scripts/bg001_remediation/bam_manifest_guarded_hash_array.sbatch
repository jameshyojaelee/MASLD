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
umask 077
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"

# Keep the source-recount hash-array producer byte-identical.  This frozen outer
# wrapper supplies the explicit private-file umask, then delegates every
# scientific/provenance operation to that exact producer.
exec bash "$SCRIPT_DIR/bam_manifest_hash_array.sbatch"
