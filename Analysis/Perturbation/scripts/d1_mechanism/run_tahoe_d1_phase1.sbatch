#!/bin/bash
#SBATCH --job-name=Tahoe
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d1_mechanism/logs/tahoe_d1_phase1_%j.log

# D1 Mechanism Phase 1: Tahoe-100M pseudobulk DE lookup.
# For each druggable hit, map gene → drug → pseudobulk DE (subset of 1026 shards)
# → emit top-100 affected downstream genes by |log2FC|.
# Set TAHOE_N_SHARDS to override default (50 ≈ 5% coverage). Use "ALL" for full.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_tahoe

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export TAHOE_N_SHARDS="${TAHOE_N_SHARDS:-50}"

set -u

python Analysis/Perturbation/scripts/d1_mechanism/tahoe_runner.py \
  --modality zero_shot \
  --context all

echo "[D1/Tahoe] runner finished"
