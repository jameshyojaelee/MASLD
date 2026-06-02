#!/bin/bash
#SBATCH --job-name=liver_zenodo_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/zenodo_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/zenodo_coloc_%j.err

set -euo pipefail
echo "=== Zenodo COLOC Integration ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/37_integrate_zenodo_coloc.R

echo "=== Complete ==="
echo "End: $(date)"
