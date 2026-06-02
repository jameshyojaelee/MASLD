#!/bin/bash
#SBATCH --job-name=st-corpus-prep
#SBATCH --partition=io
#SBATCH --time=48:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/state_train_v1/logs/corpus_prep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/state_train_v1/logs/corpus_prep_%j.err

# Build the cell-load-compatible ST training corpus from GEARS-processed
# Replogle K562 essential + Norman 2019 h5ads. CPU-only, IO-heavy.
set -eo pipefail   # `set -u` deferred — micromamba activate scripts reference
                   # unbound vars (MKL_INTERFACE_LAYER) and trip nounset early.

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$PROJECT_ROOT/.envs/perturbation_state"

# Strict env-isolation per project policy
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
set -u   # safe to enable now that env is active

cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/shared/env_setup"

python assemble_st_train_corpus_v1.py \
    --datasets replogle_k562_essential norman

echo "[corpus-prep] complete"
ls -lh "$PROJECT_ROOT/data/perturbation/state_train_v1/datasets/"*/
