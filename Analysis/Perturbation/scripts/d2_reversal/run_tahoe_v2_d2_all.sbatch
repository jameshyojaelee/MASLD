#!/bin/bash
#SBATCH --job-name=Tahoe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Perturbation/scripts/d2_reversal/logs/Tahoe-v2_%j.log

# Tahoe v2 D2 reversal — full 1885-gene invocab subset, --context all.
# Replaces the broken --wrap-with-shebang attempts (slurm wraps run under
# /bin/sh which rejects `set -o pipefail`).

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true
export CUDA_VISIBLE_DEVICES=""
eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_state
set -u

python Analysis/Perturbation/scripts/d2_reversal/tahoe_reversal_v2_runner.py \
  --modality zero_shot --context all --reference ref_a \
  --hits Analysis/Perturbation/data/hits/d2_reversal_hits_invocab.csv

echo "[D2/Tahoe v2] runner finished"
