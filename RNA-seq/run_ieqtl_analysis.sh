#!/bin/bash
#SBATCH --job-name=liver_ieqtl_analysis
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ieqtl_analysis_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ieqtl_analysis_%j.err

set -euo pipefail
echo "=== ieQTL Analysis ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/38_ieqtl_analysis.R

echo "=== Complete ==="
echo "End: $(date)"
