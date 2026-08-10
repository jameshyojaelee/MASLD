#!/usr/bin/env bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 2; }
SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
MANIFEST="$BG001_RUN_ROOT/manifests/bam_manifest.tsv"
mkdir -p "$BG001_RUN_ROOT/logs"
exec > >(tee "$BG001_RUN_ROOT/logs/merge_validate.out") 2> >(tee "$BG001_RUN_ROOT/logs/merge_validate.err" >&2)

eval "$(micromamba shell hook --shell=bash)"
set +u
micromamba activate rnaseq
set -u
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT" --require-bam-frozen
[[ -f "$BG001_RUN_ROOT/BASELINE_FROZEN" ]] || { echo "Baseline is not frozen" >&2; exit 64; }
Rscript --vanilla "$SCRIPT_DIR/merge_gse213621.R" "$BG001_RUN_ROOT" "$MANIFEST"
STAGE="$BG001_RUN_ROOT/counts/GSE213621/.merge_staging"
python3 "$SCRIPT_DIR/validate_featurecounts.py" \
  --dataset GSE213621 \
  --manifest "$MANIFEST" \
  --counts "$STAGE/gene_counts.txt" \
  --summary "$STAGE/gene_counts.txt.summary" \
  --log "$STAGE/featureCounts.log" \
  --environment "$STAGE/environment.txt" \
  --expected-samples 367 \
  --output "$STAGE/validation.json"
python3 "$SCRIPT_DIR/publish_gse213621.py" --run-root "$BG001_RUN_ROOT"

python3 "$SCRIPT_DIR/finalize_recount.py" --run-root "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/make_count_overrides.py" \
  --run-root "$BG001_RUN_ROOT" \
  --count-manifest "$BG001_RUN_ROOT/manifests/count_manifest.tsv" \
  --baseline-overrides "$BG001_RUN_ROOT/manifests/read_count_overrides.tsv" \
  --output "$BG001_RUN_ROOT/manifests/count_overrides.tsv"
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT" --require-bam-frozen
python3 "$SCRIPT_DIR/verify_recount_artifacts.py" --run-root "$BG001_RUN_ROOT" --seal
python3 "$SCRIPT_DIR/verify_recount_artifacts.py" --run-root "$BG001_RUN_ROOT"
