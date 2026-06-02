#!/bin/bash
#SBATCH --job-name=STATE
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d2_reversal/logs/state_d2_phase1_%j.log

# D2 Reversal Phase 1: STATE zero-shot — SE-600M embedding-space disease-axis
# proxy. Per-gene reversal score = - cosine(emb(gene) - pop_avg, disease_axis).
# CPU-only (torch 2.5.1+cu124 vs cluster GPU SM mismatch); SE-600M is 2.86 GB.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
export CUDA_VISIBLE_DEVICES=""
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_state

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

set -u

REFERENCE="${REFERENCE:-ref_a}"
MAX_HITS="${MAX_HITS:-500}"

python Analysis/Perturbation/scripts/d2_reversal/state_reversal_runner.py \
  --modality zero_shot \
  --context all \
  --reference "$REFERENCE" \
  --max-hits "$MAX_HITS"

echo "[D2/STATE] runner finished"
