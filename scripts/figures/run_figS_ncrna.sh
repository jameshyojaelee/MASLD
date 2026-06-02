#!/bin/bash
#SBATCH --job-name=figS_ncrna
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/figS_ncrna_%j.out
#SBATCH --error=logs/figS_ncrna_%j.err

echo "=== figS_ncrna ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p scripts/figures/logs

Rscript scripts/figures/figS_ncrna.R

echo "Exit code: $?"
echo "End: $(date)"
