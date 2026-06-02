#!/bin/bash
#SBATCH --job-name=pseudobulk_de
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/pseudobulk_de_%j.log

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="${BASE}"

mkdir -p logs

echo "=== Pseudobulk DE Analysis ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript "${BASE}/Analysis/SingleCell/scripts/pseudobulk_de.R"

echo "=== Done: $(date) ==="
