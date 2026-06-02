#!/bin/bash
#SBATCH --job-name=ctwas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/ctwas_%j.out
#SBATCH --error=RNA-seq/logs/ctwas_%j.err

# Phase 2B-i: cTWAS with GTEx Liver weights
# Submit: GWAS_NAME=ghodsian sbatch RNA-seq/run_ctwas.sbatch

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

GWAS_NAME="${GWAS_NAME:-ghodsian}" Rscript RNA-seq/64_ctwas_gtex.R
