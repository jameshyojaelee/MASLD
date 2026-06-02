#!/bin/bash
#SBATCH --job-name=B1_zonation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/B1_zonation_%j.out
#SBATCH --error=logs/B1_zonation_%j.err

echo "=== Module B1: Hepatocyte Zonation Classification ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 43_zonation_classification.R

echo "End: $(date)"
echo "Exit code: $?"
