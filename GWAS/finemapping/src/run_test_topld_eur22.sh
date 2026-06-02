#!/bin/bash -l
#SBATCH --job-name=test_topld_eur22
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=4
#SBATCH --mem=220G
#SBATCH --time=04:00:00
#SBATCH --output=GWAS/finemapping/logs/test_topld_eur22_%j.out
#SBATCH --error=GWAS/finemapping/logs/test_topld_eur22_%j.err

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript GWAS/finemapping/src/build_topld_blocks.R EUR 22
