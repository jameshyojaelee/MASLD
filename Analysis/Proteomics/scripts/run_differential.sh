#!/bin/bash
#SBATCH --job-name=proteo_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/differential_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/differential_v2_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PROTEO_DIR="${BASE}/Analysis/Proteomics"

mkdir -p "${PROTEO_DIR}/logs"

echo "=== Differential Proteomics v2 (corrected contrasts) ==="
echo "Time: $(date)"
echo "Node: $(hostname)"

# Activate environment (guard against `set -u`: conda activation scripts reference
# unbound vars like ADDR2LINE which abort under `set -euo pipefail`).
eval "$(micromamba shell hook -s bash)"
set +u
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT="${BASE}"

Rscript "${PROTEO_DIR}/scripts/differential_proteomics.R"

echo "Done: $(date)"
