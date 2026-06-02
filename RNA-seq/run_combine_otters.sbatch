#!/bin/bash
#SBATCH --job-name=combine_otters
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/combine_otters_%j.out
#SBATCH --error=RNA-seq/logs/combine_otters_%j.err

# Combine per-chromosome OTTERS TWAS results with ACAT
set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Combining OTTERS TWAS Results ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"

python RNA-seq/63c_combine_otters_results.py

echo "=== Done: $(date) ==="
