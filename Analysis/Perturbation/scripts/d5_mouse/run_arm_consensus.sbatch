#!/bin/bash
#SBATCH --job-name=D5-consensus
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/d5_consensus_%j.out
#SBATCH --error=logs/d5_consensus_%j.err

# D5 Mouse cross-species consensus + concordance_scorer + backtester + atlas_integrator.

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/d5_mouse"

eval "$(micromamba shell hook --shell bash)"
micromamba activate perturbation_state
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true
set -u

python consensus.py
python concordance_scorer.py
python backtester.py
python atlas_integrator.py
echo "[D5] consensus pipeline complete"
