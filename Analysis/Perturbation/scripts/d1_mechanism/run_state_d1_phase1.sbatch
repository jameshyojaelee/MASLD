#!/bin/bash
#SBATCH --job-name=STATE
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=4:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d1_mechanism/logs/state_d1_phase1_%j.log

# D1 Mechanism Phase 1: STATE zero-shot — SE-600M embedding-neighbor proxy
# for downstream gene prediction. Top-100 nearest neighbours in 2048-d embedding
# space emitted as predicted downstream-affected genes.
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

python Analysis/Perturbation/scripts/d1_mechanism/state_runner.py \
  --modality zero_shot \
  --context all

echo "[D1/STATE] runner finished"
