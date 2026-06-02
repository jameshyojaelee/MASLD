#!/bin/bash
#SBATCH --job-name=A3_convergence
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/A3_convergence_%j.out
#SBATCH --error=logs/A3_convergence_%j.err

echo "=== Module A3: Regulon-Drug-GWAS Convergence ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 42_regulon_drug_gwas_convergence.R

echo "End: $(date)"
echo "Exit code: $?"
