#!/bin/bash
#SBATCH --job-name=Tahoe
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d2_reversal/logs/tahoe_d2_phase1_%j.log

# D2 Reversal Phase 1: Tahoe-100M drug-pseudobulk vs reference LFC cosine.
# CPU-only — we read parquet shards + compute pandas/numpy aggregations; no
# Tahoe model inference (the dataset's pseudobulk DE is precomputed by Vevo
# and shipped with the 429 GB dataset).
# HepG2/C3A-only mode (shards 432-451) keeps the proxy tissue-relevant.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
export CUDA_VISIBLE_DEVICES=""
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
# Use perturbation_state env — it has pandas/numpy/pyarrow installed; the
# dedicated perturbation_tahoe env is currently broken (task #21).
micromamba activate -p .envs/perturbation_state

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

set -u

REFERENCE="${REFERENCE:-ref_a}"
MAX_HITS="${MAX_HITS:-500}"
export TAHOE_HEPG2_ONLY="${TAHOE_HEPG2_ONLY:-1}"

python Analysis/Perturbation/scripts/d2_reversal/tahoe_reversal_runner.py \
  --modality zero_shot \
  --context all \
  --reference "$REFERENCE" \
  --max-hits "$MAX_HITS"

echo "[D2/Tahoe] runner finished"
