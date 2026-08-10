#!/bin/bash
#SBATCH --job-name=susie_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/paltest_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/paltest_%j.err
set -eo pipefail
set +u
export PATH="$HOME/.local/bin:$PATH"; eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export COLOC_OUT_SUFFIX=_paltest
export GWAS_NAME=2021_34841290_NAFLD_EUR
export CHR_FILTER=22
Rscript 06_susie_coloc.R
