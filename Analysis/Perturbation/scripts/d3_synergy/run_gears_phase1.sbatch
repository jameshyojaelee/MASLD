#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_phase1_%j.log

# GEARS Phase 1: auto-trains on Replogle K562 + predicts top-100 D3 tier1 pairs.
# Hours-long run (Replogle download + 5-epoch train).

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

python Analysis/Perturbation/scripts/d3_synergy/gears_runner.py \
  --modality zero_shot --tier tier1_pairs --max-hits 100
