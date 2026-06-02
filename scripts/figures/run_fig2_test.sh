#!/bin/bash
#SBATCH --job-name=fig2_test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/fig2_test_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/fig2_test_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="${BASE}"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org")'
echo "--- Fig 1 ---"
Rscript "${BASE}/scripts/figures/fig1_compact.R"
