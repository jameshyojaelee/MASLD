#!/bin/bash
#SBATCH --job-name=alpha3_mofa
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/503_mofa_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/503_mofa_%j.err
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp
source "$HOME/.bashrc"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# Step 1: build pseudobulk TSV views (anndata available in rapids_singlecell)
micromamba activate rapids_singlecell
python -u 503_mofa_pseudobulk_prep.py
micromamba deactivate

# Step 2: fit MOFA+ per view set in R (celltype_bio has MOFA2)
micromamba activate celltype_bio
Rscript 503_mofa_pseudobulk.R
