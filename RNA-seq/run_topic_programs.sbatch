#!/bin/bash
#SBATCH --job-name=topic_programs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/topic_programs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/topic_programs_%j.err

set -euo pipefail

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs

echo "=== Topic Programs (NMF) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

Rscript "${MASLD_PROJECT_ROOT}/RNA-seq/230_topic_programs.R"

echo ""
echo "End: $(date)"
