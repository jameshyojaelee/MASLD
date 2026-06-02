#!/bin/bash
#SBATCH --job-name=sex_figs_final
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_figs_final_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/sex_figs_final_%j.err

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="${BASE}/scripts/figures"
mkdir -p "${BASE}/figures/logs"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Regenerating sex figures (post-aggregation) ==="
echo "Date: $(date)"

echo "--- figS_sex_divergent.R ---"
Rscript "${SCRIPT_DIR}/figS_sex_divergent.R"

echo "--- figS_sex_subsampling.R ---"
Rscript "${SCRIPT_DIR}/figS_sex_subsampling.R"

echo "--- fig2_disease_progression.R ---"
Rscript "${SCRIPT_DIR}/fig2_disease_progression.R"

echo "=== Done: $(date) ==="
