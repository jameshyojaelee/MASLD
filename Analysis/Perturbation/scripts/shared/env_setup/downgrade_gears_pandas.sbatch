#!/bin/bash
#SBATCH --job-name=GEARS-fix
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/gears_pandas_%j.log

# Downgrade pandas in gears env to <2 — GEARS uses Series.nonzero() which was
# removed in pandas 2.x. Conda install (NOT pip) to avoid user-site leak.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_gears

eval "$(micromamba shell hook --shell bash)"

echo "[+] Downgrading pandas + scipy in gears env"
micromamba install -p "$ENV_PATH" -c conda-forge -y 'pandas=1.5.3' 'scipy<1.13'

echo ""
echo "[+] Verify"
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -c "
import pandas, scipy, numpy
print('pandas:', pandas.__version__)
print('scipy:', scipy.__version__)
print('numpy:', numpy.__version__)
import gears; print('gears OK')
from gears import GEARS, PertData; print('GEARS+PertData OK')
"
