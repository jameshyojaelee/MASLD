#!/bin/bash
#SBATCH --job-name=B5_biomart
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs/B5_biomart_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs/B5_biomart_%j.err

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs

source ~/.bashrc
micromamba activate rnaseq
set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
Rscript scripts/B5_biomart_pin_ensembl112.R
