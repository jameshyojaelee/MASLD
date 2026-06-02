#!/bin/bash
#SBATCH --job-name=figS15_model_comparison
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS15_model_comparison_%j.out
#SBATCH --error=logs/figS15_model_comparison_%j.err

set -euo pipefail

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript scripts/figures/assemble_S15_model_comparison.R

echo "Done: $(date)"
