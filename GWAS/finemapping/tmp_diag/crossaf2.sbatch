#!/bin/bash
#SBATCH --job-name=susieR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/crossaf2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/crossaf2_%j.err
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping
set -u
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/crossstudy_af2.R
