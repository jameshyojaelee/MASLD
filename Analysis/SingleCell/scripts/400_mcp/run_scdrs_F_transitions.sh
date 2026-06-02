#!/bin/bash
#SBATCH --job-name=scdrs_Ftrans
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/scdrs_Ftrans.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/scdrs_Ftrans.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate spatial

echo "[scdrs] host=$(hostname)  cwd=$(pwd)"
echo "[scdrs] python=$(which python)"
python --version
python -c "import scdrs; print('scdrs', scdrs.__version__)"

python run_scdrs_F_transitions.py

echo "[scdrs] DONE."
