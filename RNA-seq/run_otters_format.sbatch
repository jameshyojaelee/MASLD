#!/bin/bash
#SBATCH --job-name=otters_format
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/otters_format_%j.out
#SBATCH --error=RNA-seq/logs/otters_format_%j.err

set -euo pipefail

echo "=== OTTERS Format Broadaway eQTL ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Create logs directory if needed
mkdir -p logs

# Activate environment
# NOTE: rnaseq env has a pandas/numpy C-extension conflict (pandas 2.3.x).
# snapatac2 has the same versions but with compatible builds.
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2

# Project root
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

# Run the formatting script
python "${PROJECT_ROOT}/RNA-seq/61_otters_format_broadaway.py" \
    --project-root "${PROJECT_ROOT}" \
    --min-n 100 \
    --chromosomes 1-22

echo ""
echo "=== Finished: $(date) ==="
