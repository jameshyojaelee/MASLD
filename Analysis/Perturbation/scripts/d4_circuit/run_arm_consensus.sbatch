#!/bin/bash
#SBATCH --job-name=D4-consensus
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/d4_consensus_%j.out
#SBATCH --error=logs/d4_consensus_%j.err

# D4 Circuit consensus + backtester + atlas_integrator post-pass.
# Also runs lr_curator as a pre-step (was racing in old run_arm.sbatch P1-D4-5).

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/d4_circuit"

eval "$(micromamba shell hook --shell bash)"
micromamba activate perturbation_ccc
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true
set -u

# LR-pair curation (idempotent; cheap)
python lr_curator.py || echo "[D4] lr_curator failed but continuing"

python consensus.py
python backtester.py
python atlas_integrator.py
echo "[D4] consensus pipeline complete"
