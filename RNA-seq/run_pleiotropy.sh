#!/bin/bash
#SBATCH --job-name=masld_pleiotropy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=logs/pleiotropy_%j.out
#SBATCH --error=logs/pleiotropy_%j.err

set -euo pipefail
echo "=== Pleiotropy Analysis ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs
Rscript RNA-seq/39_pleiotropy_analysis.R

echo "=== Complete ==="
echo "End: $(date)"
