#!/bin/bash
#SBATCH --job-name=test_rsc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/test_rsc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/test_rsc_%j.err

set -euo pipefail

echo "=== rapids_singlecell Validation ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
echo ""

# Show GPU info
nvidia-smi

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

echo ""
echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Run validation
python /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/test_rapids_singlecell.py

echo ""
echo "Done: $(date)"
