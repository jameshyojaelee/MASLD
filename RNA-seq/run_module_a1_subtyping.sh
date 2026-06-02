#!/bin/bash
#SBATCH --job-name=A1_subtyping
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/A1_subtyping_%j.out
#SBATCH --error=logs/A1_subtyping_%j.err

echo "=== Module A1: Molecular Subtyping ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 44_molecular_subtyping.R

echo "End: $(date)"
echo "Exit code: $?"
