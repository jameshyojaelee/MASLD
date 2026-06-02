#!/bin/bash
#SBATCH --job-name=cnmf_seed2024_prep
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_seed2024_prep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_seed2024_prep_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

python 420_run_cnmf.py prepare \
  --input atlas_cnmf_global.h5ad \
  --name global_seed2024 \
  --k 16 \
  --n-iter 200 \
  --num-hvg 2500 \
  --seed 2024

echo "[420.prepare seed2024] DONE."
