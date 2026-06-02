#!/bin/bash
#SBATCH --job-name=B2_ferroptosis
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/B2_ferroptosis_%j.out
#SBATCH --error=logs/B2_ferroptosis_%j.err

echo "=== Module B2: Ferroptosis/Lipotoxicity Program ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 41_ferroptosis_lipotoxicity.R

echo "End: $(date)"
echo "Exit code: $?"
