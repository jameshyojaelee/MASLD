#!/bin/bash
#SBATCH --job-name=B5_nmf_seeds
#SBATCH --partition=bigmem
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs/B5_nmf_seeds_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs/B5_nmf_seeds_%j.err

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs

source ~/.bashrc
micromamba activate rnaseq
set -eo pipefail

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export NMF_RUNS=50

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
Rscript scripts/B5_nmf_seed_robustness.R
