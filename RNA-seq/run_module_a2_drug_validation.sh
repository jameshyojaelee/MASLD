#!/bin/bash
#SBATCH --job-name=A2_drug_validation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/A2_drug_validation_%j.out
#SBATCH --error=logs/A2_drug_validation_%j.err

echo "=== Module A2: Drug-Target Clinical Validation ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 40_drug_target_validation.R

echo "End: $(date)"
echo "Exit code: $?"
