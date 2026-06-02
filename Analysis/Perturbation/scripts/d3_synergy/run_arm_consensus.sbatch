#!/bin/bash
#SBATCH --job-name=D3-consensus
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/d3_consensus_%j.out
#SBATCH --error=logs/d3_consensus_%j.err

# D3 Synergy consensus + backtester + atlas_integrator post-pass.

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/d3_synergy"

eval "$(micromamba shell hook --shell bash)"
micromamba activate perturbation_state
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true
set -u

python consensus.py
python backtester.py
python atlas_integrator.py
echo "[D3] consensus pipeline complete"
