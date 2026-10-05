#!/usr/bin/env bash
#SBATCH --job-name=1593
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1593_verify_lineage_losses_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1593_verify_lineage_losses_%j.err
set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
cd "$TASK_ROOT"
: "${LINEAGE_VERIFY_MODEL_SUMMARY_SHA256:?set the completed model summary SHA256}"
: "${LINEAGE_VERIFY_MODEL_ARTIFACTS_SHA256:?set the completed model artifacts SHA256}"
exec /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse296875_verify_lineage_losses.py" \
  --expected-script-sha256 ab002afce91c85177cb25ffda114084aa1fba09300d37f8cba7c867885109219 \
  --expected-summary-sha256 "$LINEAGE_VERIFY_MODEL_SUMMARY_SHA256" \
  --expected-artifacts-sha256 "$LINEAGE_VERIFY_MODEL_ARTIFACTS_SHA256"
