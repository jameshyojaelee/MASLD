#!/bin/bash
#SBATCH --job-name=gen_gene_profiles
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gen_gene_profiles_%j.log
#SBATCH --error=logs/gen_gene_profiles_%j.log

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/masld-atlas-v2
mkdir -p scripts/logs

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== generate_gene_profiles.py ==="
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Start: $(date)"

python scripts/generate_gene_profiles.py --output-dir public/data

echo "End: $(date)"
echo "=== done ==="
