#!/bin/bash
#SBATCH --job-name=cnmf_seed2024_cons
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_seed2024_cons_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_seed2024_cons_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

python 420_run_cnmf.py consensus \
  --name global_seed2024 \
  --k 16 \
  --density-threshold 0.03

echo "[420.consensus seed2024] DONE."
