#!/bin/bash
#SBATCH --job-name=broadaway
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/logs/63_broadaway_bench_prep_%j.log

set -eo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$BASE/GWAS/finemapping/src/logs"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd "$BASE/GWAS/finemapping/src"
Rscript 63_broadaway_bench_prep.R
