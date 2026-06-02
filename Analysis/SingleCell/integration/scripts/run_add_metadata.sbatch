#!/bin/bash
#SBATCH --job-name=add_sc_metadata
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=Analysis/SingleCell/integration/logs/add_metadata_%j.log

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Add Sample Metadata to SC Atlas ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

python3 Analysis/SingleCell/integration/scripts/add_sample_metadata.py

echo "=== Done: $(date) ==="
