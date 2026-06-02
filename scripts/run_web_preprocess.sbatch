#!/bin/bash
#SBATCH --job-name=web_preprocess
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/web_preprocess_%j.log

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial
export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Web preprocessing: $(date) ==="
python masld-atlas-v2/scripts/preprocess_atlas_data.py
echo "=== Completed: $(date) ==="
