#!/bin/bash
#SBATCH --job-name=S5_fstage_ext
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/S5/fstage_ext_%j.out
#SBATCH --error=logs/S5/fstage_ext_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p logs/S5
source ~/.bashrc
micromamba activate rapids_singlecell
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export S5_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
python Analysis/SingleCell/scripts/354_fstage_external_check.py
