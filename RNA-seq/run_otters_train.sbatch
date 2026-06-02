#!/bin/bash
#SBATCH --job-name=otters_train
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/otters_train_chr%a_%j.out
#SBATCH --error=RNA-seq/logs/otters_train_chr%a_%j.err
#SBATCH --array=1-22

# Phase 2A Step 2: Train OTTERS weights per chromosome
# Submit: sbatch RNA-seq/run_otters_train.sbatch

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

# Load PLINK2 for pgen→bed conversion
module load PLINK/2.0a5.13

# Activate Python env with OTTERS deps
eval "$(micromamba shell hook -s bash)"
micromamba activate causal_inference

# Prevent user-level packages from interfering with conda env
export PYTHONNOUSERSITE=1

bash RNA-seq/62_otters_train_weights.sh
