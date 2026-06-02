#!/bin/bash
#SBATCH --job-name=allele_conc_56c
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=../logs/56c_allele_conc_%j.out
#SBATCH --error=../logs/56c_allele_conc_%j.err
set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate motifbreakr

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
echo "Running 56c_allele_concordance.R on $(hostname) at $(date)"
Rscript src/56c_allele_concordance.R
echo "Done at $(date)"
