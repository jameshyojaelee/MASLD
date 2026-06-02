#!/bin/bash
#SBATCH --job-name=arc-state-nd
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=30:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/install_arc_state_nodeps_%j.log

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_state

set -u

# Use --no-deps to prevent torch / numpy / etc downgrade
echo "[+] pip install arc-state --no-deps"
pip install --no-cache-dir --no-deps arc-state 2>&1 | tail -3

# Check what extras arc-state needs and install them via conda (avoids torch chaos)
echo ""
echo "[+] Identify missing transitive deps"
python -c "
try:
    import state
    print('OK state imported')
except ImportError as e:
    print(f'IMPORT FAILED: {e}')
" 2>&1 | tail -3

echo ""
echo "[+] Verify torch still at 2.5.1"
python -c "import torch; print('torch:', torch.__version__, 'cuda:', torch.version.cuda)"

echo "[done]"
