#!/bin/bash
#SBATCH --job-name=seekr_kmer
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/seekr_%j.out
#SBATCH --error=logs/seekr_%j.err

# L6 Seekr k-mer profile layer for cross-species lncRNA orthology
# Runs Seekr at k=5 and k=6 on gene-level representatives (~33k mouse x ~36k human)

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/seekr_env

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/scripts/ortholog_pipeline

echo "=== Seekr L6 k-mer orthology pipeline ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-8}"
echo "Memory: $(free -h | grep Mem | awk '{print $2}')"
echo ""

python run_seekr.py --k-values 5 6 --top-n 5 --batch-size 5000

echo ""
echo "=== Done ==="
echo "Date: $(date)"
