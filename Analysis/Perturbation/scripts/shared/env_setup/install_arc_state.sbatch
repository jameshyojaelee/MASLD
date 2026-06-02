#!/bin/bash
#SBATCH --job-name=arc-state
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/install_arc_state_%j.log

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_state

set -u

echo "=== Try pip install arc-state ==="
pip install --no-cache-dir arc-state 2>&1 | tail -8 || echo "PyPI arc-state failed"

echo ""
echo "=== Verify ==="
python -c "
try:
    import state
    print('state:', dir(state)[:10])
except Exception as e: print('state FAIL:', e)
try:
    from state import State, StateEmbeddingModel
    print('OK State + StateEmbeddingModel imported')
except Exception as e: print('State import FAIL:', str(e)[:150])
try:
    import arc_state
    print('arc_state OK:', arc_state.__version__)
except Exception as e: print('arc_state FAIL:', str(e)[:100])
"

echo "[done]"
