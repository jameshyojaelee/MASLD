#!/bin/bash
#SBATCH --job-name=disease_sig_score
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=logs/disease_sig_scoring_%j.out
#SBATCH --error=logs/disease_sig_scoring_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell
mkdir -p logs

echo "=== Disease Signature Scoring Pipeline ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-16}"
echo "Memory: ${SLURM_MEM_PER_NODE:-200G}"

micromamba run -n rapids_singlecell python scripts/disease_signature_scoring.py

echo "=== Done: $(date) ==="
echo "Output files:"
ls -lh results_gpu_v2/disease_signatures/ 2>/dev/null || echo "No output directory found"
echo "Figures:"
ls -lh results_gpu_v2/disease_signatures/figures/ 2>/dev/null || echo "No figures directory found"
