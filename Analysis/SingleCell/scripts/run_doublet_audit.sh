#!/bin/bash
#SBATCH --job-name=S5_dbl_audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/S5/doublet_audit_%j.out
#SBATCH --error=logs/S5/doublet_audit_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p logs/S5
source ~/.bashrc
micromamba activate rnaseq
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export S5_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
python Analysis/SingleCell/scripts/355_doublet_audit.py
