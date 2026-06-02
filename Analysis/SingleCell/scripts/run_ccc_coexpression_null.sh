#!/bin/bash
#SBATCH --job-name=S5_ccc_null
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/S5/ccc_null_%j.out
#SBATCH --error=logs/S5/ccc_null_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p logs/S5
source ~/.bashrc
micromamba activate rapids_singlecell
# Read inputs from main project; write outputs to worktree (via explicit paths)
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export S5_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
export CCC_N_DONORS=${CCC_N_DONORS:-30}
export CCC_N_NULLS=${CCC_N_NULLS:-10}
python Analysis/SingleCell/scripts/353_ccc_coexpression_null.py
