#!/bin/bash
#SBATCH --job-name=D1-consensus
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/d1_consensus_%j.out
#SBATCH --error=logs/d1_consensus_%j.err

# D1 Mechanism arm — consensus + backtester + atlas_integrator post-pass.
# Submitted via run_arm_driver.sh with --dependency=afterany:$ARRAY_JOB
# (P0-X1 fix: runs regardless of per-runner success/failure).

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
ARM_DIR="$PROJECT_ROOT/Analysis/Perturbation/scripts/d1_mechanism"
cd "$ARM_DIR"

# Consensus runs in the STATE env (has pydantic + pandas; lightest of the gpu envs)
eval "$(micromamba shell hook --shell bash)"
micromamba activate perturbation_state

export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true

set -u

echo "[D1/consensus] running consensus.py"
python consensus.py

echo "[D1/backtest] running backtester.py"
python backtester.py

echo "[D1/atlas] running atlas_integrator.py"
python atlas_integrator.py

echo "[D1] consensus pipeline complete"
