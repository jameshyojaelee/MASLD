#!/bin/bash
#SBATCH --job-name=STATE
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=4:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d5_mouse/logs/state_mouse_phase1_%j.log

# D5 Mouse Phase 1: STATE-mouse zero-shot — SE-600M cross-species embedding
# concordance on the Conserved_Core 1,180 mouse-human ortholog rows.
# CPU-only (torch 2.5.1+cu124 vs cluster GPU SM mismatch); SE-600M is 2.86 GB.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
export CUDA_VISIBLE_DEVICES=""
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_state

# After activation, $CONDA_PREFIX is set to the env path. Wire LD path so the
# env's own libstdc++ / CUDA-stub libs win over system libs.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

set -u

python Analysis/Perturbation/scripts/d5_mouse/state_mouse_runner.py \
  --modality zero_shot \
  --max-hits 1180

echo "[D5/STATE] runner finished"
