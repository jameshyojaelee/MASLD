#!/bin/bash
#SBATCH --job-name=figS_batch_overlap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/figS_batch_overlap_%j.out
#SBATCH --error=scripts/figures/logs/figS_batch_overlap_%j.err

set -eo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

set +u
source ~/.bashrc
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT="${BASE}"

echo "[$(date '+%H:%M:%S')] Running figS_batch_correction_overlap.R"
Rscript scripts/figures/figS_batch_correction_overlap.R
echo "[$(date '+%H:%M:%S')] DONE"
