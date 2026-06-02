#!/bin/bash
#SBATCH --job-name=otters_twas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/otters_twas_%j.out
#SBATCH --error=RNA-seq/logs/otters_twas_%j.err

# Phase 2A Step 3: Run OTTERS TWAS per GWAS
# Submit: GWAS_NAME=ghodsian sbatch RNA-seq/run_otters_twas.sbatch

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

module load PLINK/2.0a5.13

eval "$(micromamba shell hook -s bash)"
micromamba activate causal_inference

# Prevent user-level packages from interfering with conda env
export PYTHONNOUSERSITE=1

bash RNA-seq/63_otters_twas.sh
